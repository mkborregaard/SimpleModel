# The model layer. Everything is expressed in terms of a single SpatialEcology
# `Assemblage` that is raster-backed (its `Locations` wrap a `RasterData`):
#
#   * the species-by-sites occurrence matrix replaces the old `Species.ranges`,
#   * the domain mask and per-site cell indices replace `Environment.mask`/`.inds`,
#   * the two environmental PCA axes live as the `:pca1` / `:pca2` site statistics.
#
# So the `Environment` and `Species` structs are gone — the Assemblage *is* the
# data object, and most of the old bookkeeping is provided by SpatialEcology.

using Rasters
using Statistics
using StatsBase: fit, Histogram
using GeoInterface; const GI = GeoInterface
using LibGEOS
import ConcaveHull          # `import`, not `using`: ConcaveHull also exports `area`,
                            # which would clash with LibGEOS's
using SpatialEcology
import SpatialEcology: places, getcoords
import GeometryBasics      # `import`, not `using`: GeometryBasics also exports
                           # `coordinates`, which would clash with SpatialEcology's

# ---------------------------------------------------------------------------
# Bridge accessors: the raster-backed ingredients the model needs from an
# Assemblage. These are the only places that reach into the Locations.
# ---------------------------------------------------------------------------

"The Bool domain raster the assemblage is defined on."
domain(asm) = getcoords(places(asm)).mask

"One `CartesianIndex` into the domain raster per site (in site order)."
cellindices(asm) = getcoords(places(asm)).cellinds

"A site statistic (e.g. `:pca1`) scattered back onto the domain as a Raster."
pca_map(asm, var::Symbol) = to_raster(asm[!, var], asm)

"A single species' range as a Bool raster."
species_raster(asm, sp) = to_raster(Vector(getspecies(asm, sp)), asm)

"The climate-space (PCA) coordinates of the sites where species `sp` occurs."
function climate(asm, sp)
    s = occupied(asm, sp)
    asm[!, :pca1][s], asm[!, :pca2][s]
end

# ---------------------------------------------------------------------------
# Climate-space geometry helpers
# ---------------------------------------------------------------------------

function points_to_geo(xs, ys)
    length(xs) < 1 && return GI.MultiPoint([(0, 0)])
    GI.MultiPoint(collect(zip(xs, ys)))
end
points_to_geo(points) = GI.MultiPoint(points)

# Bounding box of all sites in PCA space, as ((xmin, xmax), (ymin, ymax)).
climate_bbox(asm) = (extrema(asm[!, :pca1]), extrema(asm[!, :pca2]))

# Concave hull of all sites in PCA space, as a GeoInterface polygon.
function climate_hull(asm; k = 40)
    pca1, pca2 = asm[!, :pca1], asm[!, :pca2]
    cc = ConcaveHull.concave_hull(ConcaveHull.KDTree(vcat(pca1', pca2'), reorder = false), k)
    GI.Polygon([GI.LinearRing([cc.vertices; [first(cc.vertices)]])])
end

# ---------------------------------------------------------------------------
# Per-species summaries in climate / geographic space
# ---------------------------------------------------------------------------

# Centroid in climate space of a species' occurrences.
function get_centroid(asm, sp)
    x, y = climate(asm, sp)
    mean(x), mean(y)
end

# Geographic centroids (x, y) of every species, in species order.
function geocentroids(asm)
    co = coordinates(asm)
    map(1:nspecies(asm)) do i
        s = occupied(asm, i)
        isempty(s) ? (NaN, NaN) : (mean(@view co[s, 1]), mean(@view co[s, 2]))
    end
end

# Shape descriptors of each species' climate-space occurrences: the correlation
# of the two axes and the spread along each.
function find_range_shapes(asm)
    cors = Float64[]; xrange = Float64[]; yrange = Float64[]
    for sp in speciesnames(asm)
        x, y = climate(asm, sp)
        if isempty(x)
            push!(cors, 0); push!(xrange, 0); push!(yrange, 0)
        else
            push!(cors, cor(x, y))
            push!(xrange, maximum(x) - minimum(x))
            push!(yrange, maximum(y) - minimum(y))
        end
    end
    cors, xrange, yrange
end

# Convex hull of a species' occurrences in climate space, and its area.
gethull(asm, sp) = convexhull(points_to_geo(climate(asm, sp)...))
hullarea(asm, sp) = LibGEOS.area(gethull(asm, sp))

# ---------------------------------------------------------------------------
# Ellipse fitting and sampling in climate space
# ---------------------------------------------------------------------------

# Fraction of a sampled ellipse's area that falls within `polygon` (the climate
# hull) — used to keep randomly placed ellipses inside the realised climate.
function overlap(el::Ellipse, polygon; n = 100)
    ellipse_points = points_to_geo(GeometryBasics.coordinates(el, n))
    ellipse_poly = GI.Polygon([GI.LinearRing(GI.getpoint(ellipse_points))])
    LibGEOS.area(LibGEOS.intersection(polygon, ellipse_poly)) / LibGEOS.area(ellipse_poly)
end

# Per-point inverse density in climate space: 1 / (number of occurrences sharing
# the same `binsize` grid cell). Used to down-weight oversampled climates.
function makeweights(xs, ys, binsize = 0.1)
    makebins(i) = floor(minimum(i) - binsize):binsize:ceil(maximum(i) + binsize)
    points_in_cell(x, y, hist) = hist.weights[findfirst(>(x), hist.edges[1]) - 1, findfirst(>(y), hist.edges[2]) - 1]

    hist = fit(Histogram, (xs, ys), (makebins(xs), makebins(ys)))
    [1 / points_in_cell(xs[i], ys[i], hist) for i in eachindex(xs, ys)]
end

# Fit an elliptical niche to a species' occurrences in climate space. `method`
# chooses the point-trust / quality-control strategy (see TrustMethod in
# ellipse.jl) used to decide which occurrences are trusted true presences; the
# shape is then fit to those points and the boundary sized to enclose the `cover`
# quantile of them (cover = nothing sizes by `sigma` core SDs instead). Swap
# `method` to compare quality-control techniques. `weightmap`, if given, is a
# per-site weight vector (e.g. from `makeweights`).
function fitellipse(asm, sp, sigma = 2;
        weightmap = nothing, method::TrustMethod = MCDTrim(), cover = 1.0)
    xs, ys = climate(asm, sp)
    length(xs) < 3 && return Ellipse(0, 0, 0, 0, 0)
    weight = isnothing(weightmap) ? ones(length(xs)) : weightmap[occupied(asm, sp)]
    fit(Ellipse, xs, ys, sigma; weight, method, cover)
end

# The climate-space occurrences of a species together with their per-point trust
# under a given method (1 = trusted true presence, 0 = discarded) — for
# diagnostics and plotting.
function species_trust(asm, sp; weightmap = nothing, method::TrustMethod = MCDTrim())
    xs, ys = climate(asm, sp)
    weight = isnothing(weightmap) ? ones(length(xs)) : weightmap[occupied(asm, sp)]
    xs, ys, trust(method, xs, ys; weight)
end

# Pick a centre for a random ellipse: either a real occupied grid cell (so the
# ellipse is anchored on realised climate) or a uniform point in the bbox.
function pick_ellipse_center(asm; on_real_point = false)
    pca1, pca2 = asm[!, :pca1], asm[!, :pca2]
    if on_real_point
        i = rand(1:length(pca1))
        return pca1[i], pca2[i]
    end
    bb = climate_bbox(asm)
    rescale(rand(), first(bb)...), rescale(rand(), last(bb)...)
end

# Draw a random ellipse of the given area that sits mostly inside the realised
# climate (`chull` is the assemblage's climate hull, see `climate_hull`).
function sample_ellipse(harea = 1; asm, chull, max_iter = 1000, on_real_point = false)
    el = Ellipse(0, 0, 0, 0, 0)
    harea == 0 && return el
    ovrlp = 0
    failsafe = 0
    while ovrlp < 0.8 && (failsafe += 1) < max_iter
        el = rand(Ellipse, pick_ellipse_center(asm; on_real_point)..., area = harea)
        ovrlp = overlap(el, chull)
    end
    failsafe == max_iter && @show harea
    el
end
