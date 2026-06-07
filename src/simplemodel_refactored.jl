# Draft: how simplemodel.jl functions change when dispatching on Assemblage.
# Not all functions are shown — only those whose signatures meaningfully change.
# Functions that operate purely on Ellipse, pca vectors, etc. are unchanged.

using SpatialEcology

# do_map: now takes an Assemblage instead of a bare mask
do_map(x::AbstractVector, asm::Assemblage) = to_raster(x, asm)

# get_climate: extract pca coords for sites where a species is present
function get_climate(species::String, asm::Assemblage)
    sp_row  = findfirst(==(species), speciesnames(asm))
    present = findnz(SpatialEcology.occurrences(asm)[sp_row, :])[1]
    ss = sitestats(asm)
    ss.pca1[present], ss.pca2[present]
end

# get_centroid: mean pca position of a species
function get_centroid(species::String, asm::Assemblage)
    xs, ys = get_climate(species, asm)
    mean(xs), mean(ys)
end

# allcentroids and rangesizes become one-liners using SpatialEcology's API:
#   allcentroids = [get_centroid(sp, asm) for sp in speciesnames(asm)]
#   rangesizes   = vec(SpatialEcology.occupancy(asm))   # sites per species

# find_range_shapes: unchanged logic, new dispatch
function find_range_shapes(asm::Assemblage)
    cors, xrange, yrange = Float64[], Float64[], Float64[]
    for sp in speciesnames(asm)
        x, y = get_climate(sp, asm)
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

# gethull / hullarea: unchanged logic, new dispatch
gethull(species::String, asm::Assemblage) =
    LibGEOS.convexhull(points_to_geo(get_climate(species, asm)...))
hullarea(species::String, asm::Assemblage) =
    LibGEOS.area(gethull(species, asm))

# geocentroids: unchanged logic, now uses Assemblage
function geocentroids(name::String, asm::Assemblage)
    sp_row  = findfirst(==(name), speciesnames(asm))
    present = findnz(SpatialEcology.occurrences(asm)[sp_row, :])[1]
    site_inds = inds(asm)[present]   # (row, col) tuples from sitestats
    isempty(site_inds) && return (NaN, NaN)
    mean(first, site_inds), mean(last, site_inds)
end

# makeweights: unchanged — operates on pca vectors, call as makeweights(pca1(asm), pca2(asm))

# fitellipse: new dispatch
function fitellipse(species::String, asm::Assemblage, sigma = 2; weightmap = nothing)
    xs, ys = get_climate(species, asm)
    length(xs) < 3 && return Ellipse(0, 0, 0, 0, 0)
    sp_row  = findfirst(==(species), speciesnames(asm))
    present = findnz(SpatialEcology.occurrences(asm)[sp_row, :])[1]
    w = isnothing(weightmap) ? ones(length(xs)) : weightmap[present]
    fit(Ellipse, xs, ys, sigma; weight = w)
end

# pick_ellipse_center: takes EnvSummary instead of Environment for bbox/chull
function pick_ellipse_center(asm::Assemblage, env::EnvSummary; on_real_point = false)
    if on_real_point
        i = rand(1:nsites(asm))
        ss = sitestats(asm)
        return ss.pca1[i], ss.pca2[i]
    end
    rescale(rand(), first(env.bbox)...), rescale(rand(), last(env.bbox)...)
end

function sample_ellipse(harea = 1; asm::Assemblage, env::EnvSummary,
                        max_iter = 1000, on_real_point = false)
    el = Ellipse(0, 0, 0, 0, 0)
    harea == 0 && return el
    ovrlp = 0; failsafe = 0
    while ovrlp < 0.8 && (failsafe += 1) < max_iter
        el = rand(Ellipse, pick_ellipse_center(asm, env; on_real_point)..., area = harea)
        ovrlp = overlap(el, env.chull)
    end
    failsafe == max_iter && @show harea
    el
end

# map_ellipse: now returns a Raster directly via to_raster
map_ellipse(el::Ellipse, asm::Assemblage) = begin
    ss = sitestats(asm)
    to_raster([in_ellipse((ss.pca1[i], ss.pca2[i]), el) for i in 1:nsites(asm)], asm)
end

# make_continuous_range: the raster-domain BFS is unchanged internally;
# we just extract a Raster{Bool} domain from the assemblage
function make_continuous_range(el::Ellipse, asm::Assemblage)
    domain = map_ellipse(el, asm)
    GeometryBasics.area(el) == 0 && return domain
    ss = sitestats(asm)
    i  = random_point_on_ellipse(el, ss.pca1, ss.pca2)
    pt = inds(asm)[i]
    growrange(pt, domain)
end
