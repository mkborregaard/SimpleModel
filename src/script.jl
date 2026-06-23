using Plots
using GLMakie
using Rasters
using JLD2
using Colors
using SpatialEcology
import GeometryBasics

include("ellipse.jl")       # the Ellipse type and the niche-fitting trust methods
include("simplemodel.jl")   # the model layer, expressed over a SpatialEcology Assemblage
include("spread.jl")        # turning niche ellipses into contiguous geographic ranges
include("plotting.jl")      # figures

save_figures = false

###--- The data is a single raster-backed Assemblage of South American birds,
# carrying the two environmental PCA axes as the :pca1 / :pca2 site statistics.
# Loading takes a while the first time, so we cache the processed object.

datadir = "/Users/cvg147/Dropbox/Arbejde/Data"
#datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
#datadir = "C:\\Users\\cvg147\\Dropbox\\Arbejde\\Data"
asm = try
    JLD2.load(joinpath(datadir, "processed_assemblage.jld2"))["asm"]
catch
    include("prepare_data.jl")
    prepare_data(datadir, doplots = true)
end

# handles used throughout: the two climate axes, one value per site
pca1, pca2 = asm[!, :pca1], asm[!, :pca2]

###--- Exploratory data analysis

# Count total diversity (species richness per cell)
diversity = richness(asm)
Plots.heatmap(richness_raster(asm), color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/empirical_richness.png")

# Plot the diversity in climate space
f = Figure()
a, s = GLMakie.scatter(f[1,1], collect(zip(pca1, pca2)); markersize = 1, color = diversity, colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2], s)
display(f)

# plot some random species and look at their distribution
randspecs = rand(speciesnames(asm), 5)
plot_species(randspecs, asm)

## how are ranges shaped?
cors, xrange, yrange = find_range_shapes(asm)
histogram(cors)
Plots.scatter(xrange, yrange)

###--- Patterns of range and niche size in climate space

# environmental centroids and (geographic) range sizes for all species
allcentroids = [get_centroid(asm, sp) for sp in speciesnames(asm)]
rangesizes = occupancy(asm)

# geographic range size on the environmental centroids (in climate space)
overplot_pca_space(allcentroids, rangesizes, asm)

# climatic range (convex-hull area) for all species, on the same centroids
allhulls = [hullarea(asm, sp) for sp in speciesnames(asm)]
overplot_pca_space(allcentroids, allhulls, asm)

# geographical centroids for all species, plotted on the map
gc = geocentroids(asm)
overplot_geo_space(gc, rangesizes, asm)

# divide range sizes into quantiles and map each quantile's centroids
rangequants = asquantiles(rangesizes, 4)
f = Figure()
aspect = DataAspect()
axs = [Axis(f[1,1]; aspect), Axis(f[2,1]; aspect), Axis(f[1,2]; aspect), Axis(f[2,2]; aspect)]
for i in 1:4
    inds = findall(==(i), rangequants)
    Makie.plot!(axs[i], domain(asm), colormap = :Greys)
    Makie.scatter!(axs[i], gc[inds], markersize = 2)
end
f

# Fit an ellipse to a random species in pca space
plot_species_pca(rand(speciesnames(asm)), asm, 2)

# Plot 16 random species with occurrences in pca space and fitted ellipses
p = Plots.plot([
    plot_species_pca(rand(speciesnames(asm)), asm, 2) for i in 1:16]...
, size = (1200, 1200))
save_figures && savefig(p, "figures/16 species in pca space.png")

# Control for the density of points in pca space by applying a grid.
# Map the options, e.g. binsizes from 0.1 to 0.4.
map_binsize(binsize) = Plots.heatmap(1 ./ to_raster(makeweights(pca1, pca2, binsize), asm), color = cgrad(:Spectral, rev = true), title = "binsize = $binsize")
Plots.plot([map_binsize(bs) for bs in 0.1:0.1:0.4]..., size = (800, 800))
save_figures && savefig("figures/binsizes.png")

# we conclude that we need one at 0.2 - a per-site inverse-density weight vector
weightmap = makeweights(pca1, pca2, 0.2)

p = Plots.plot([
    plot_species_pca(rand(speciesnames(asm)), asm, 1.5; weightmap) for i in 1:16]...
, size = (1200, 1200))
save_figures && savefig(p, "figures/16 species controlling for point density.png")

# Compare point-trust / quality-control techniques for one species (see TrustMethod).
# Drop in new `struct MyQC <: TrustMethod` + `trust(::MyQC, xs, ys; weight)` to add more.
trustmethods = [TrustAll(), MCDTrim(keep = 0.9), ChisqFilter(p = 0.95), DensityTrim(keep = 0.9)]
compare_trust_methods(rand(speciesnames(asm)), asm, trustmethods)
save_figures && savefig("figures/trust_methods.png")

# fit elliptical niches for all species (default method = MCDTrim(); swap to compare)
emp_ellipses = [fitellipse(asm, name, 1.5) for name in speciesnames(asm)]

# show patterns of ellipse area
ares = GeometryBasics.area.(emp_ellipses)
histogram(ares)
save_figures && savefig("figures/histogram of empirical ellipse areas.png")
Plots.scatter((el -> (el.center_x, el.center_y)).(emp_ellipses), marker_z = ares, ms = 3)
save_figures && savefig("figures/PCA centroids of empirical ellipses with area as color")

# fitted-ellipse overlap vs empirical richness, in climate space
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
el_emp_point = [count(el -> in_ellipse(pt, el), emp_ellipses) for pt in zip(pca1, pca2)]
Plots.plot(
    Plots.scatter(pca1, pca2, marker_z = el_emp_point, title = "fitted ellipse overlap"),
    Plots.scatter(pca1, pca2, marker_z = diversity, title = "empirical richness")
)
save_figures && savefig("figures/empirical_ellipse_and_empirical_pca_richness.png")

# Create random ellipses with the empirical areas, kept inside the realised climate
chull = climate_hull(asm)
rand_ellipses = [sample_ellipse(harea; asm, chull, on_real_point = true) for harea in ares]

# Show 50 random ellipses
p = Plots.scatter(pca1, pca2, mc = :grey, ms = 1, msw = 0, aspect_ratio = 1, label = "")
for el in rand(rand_ellipses, 50)
    Plots.plot!(p, el, label = "")
end
p
save_figures && savefig(p, "figures/50 random ellipses.png")

# plot the modelled and empirical richness
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
elpoint = [count(el -> in_ellipse(pt, el), rand_ellipses) for pt in zip(pca1, pca2)]
Plots.plot(
    Plots.scatter(pca1, pca2, marker_z = elpoint, title = "ellipse overlap"),
    Plots.scatter(pca1, pca2, marker_z = diversity, title = "empirical richness")
)
save_figures && savefig("figures/modelled_ellipse_and_empirical_pca_richness.png")

Plots.heatmap(to_raster(float.(el_emp_point), asm), color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness_on_empirical_ellipses.png")

Plots.heatmap(to_raster(float.(elpoint), asm), color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness_on_random_ellipses.png")

###--- Continuous ranges from the spread model

# Look at one random ellipse end to end
plot_ellipse_patches(rand(eachindex(emp_ellipses)), asm, emp_ellipses, 1; weightmap)

# range patches based on the empirical ellipses
model_ranges = RasterSeries([make_continuous_range(el, asm) for el in emp_ellipses], (; name = speciesnames(asm)))
newdiv = reduce(+, model_ranges)
Plots.heatmap(newdiv, color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness from patches in empirical ellipses.png")

# range patches based on randomly placed ellipses (with the right size)
rand_model_ranges = RasterSeries([make_continuous_range(el, asm) for el in rand_ellipses], (; name = speciesnames(asm)))
rand_div = reduce(+, rand_model_ranges)
Plots.heatmap(rand_div, color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness from patches in random ellipses.png")

# Find and plot richness at the 1 degree lat/long scale, via SpatialEcology's coarsen
model_coarse = richness_raster(coarsen(Assemblage(model_ranges, domain(asm)), 6))
emp_coarse = richness_raster(coarsen(asm, 6))

Plots.default(fillcolor = cgrad(:Spectral, rev = true))
Plots.plot(
    Plots.heatmap(model_coarse, title = "modelled 1 degree richness"),
    Plots.heatmap(emp_coarse, title = "empirical 1 degree richness"), size = (900, 500)
)
save_figures && savefig("figures/coarse richness.png")

# compare the range sizes of empirical ranges and those from ellipses
rand_emp_rangesize = vec(count.(model_ranges))
Plots.scatter(rangesizes, rand_emp_rangesize)
Plots.plot!(identity, 0, 6e4, color = :black, lw = 2)

# what's the relationship between ellipse size and actual range?
emp_els_area = GeometryBasics.area.(emp_ellipses)
Plots.scatter(rangesizes, emp_els_area)

Plots.default()

# How much of its niche does a species actually occupy?
function occup_in_ellipse(el::Ellipse, asm, sp)
    s = Vector(getspecies(asm, sp))
    ins, outs = 0, 0
    for i in eachindex(pca1)
        if in_ellipse((pca1[i], pca2[i]), el)
            s[i] ? (ins += 1) : (outs += 1)
        end
    end
    ins, outs
end

allins, allouts = Int[], Int[]
for ind in eachindex(emp_ellipses)
    ins, outs = occup_in_ellipse(emp_ellipses[ind], asm, ind)
    push!(allins, ins)
    push!(allouts, outs)
end

Plots.scatter(log10.(GeometryBasics.area.(emp_ellipses) .+ 0.1), allins ./ (allins .+ allouts))
Plots.scatter(log10.(rangesizes .+ 1), log10.(allins ./ (allins .+ allouts)))
# it appears that rangesize is largely determined by how much of your niche you occupy
# or - niches are consistently overestimated (possibly more likely)
# some of the small-ranged species really occur in lots of regions - why is that?

Plots.scatter(log10.(GeometryBasics.area.(emp_ellipses) .+ 0.1), log10.(rangesizes))
# larger ranges have larger ellipses but not really that strong - I believe small ellipses are exaggerated

wd = findall(>(100), allins ./ allouts)
rangesizes[wd] # the second, number 103, is almost completely occupied. Let me take a look
plot_ellipse_patches(103, asm, emp_ellipses; weightmap)
plot_ellipse_patches(3606, asm, emp_ellipses; weightmap) # The Amazonian species tend to fully occupy their ellipse

###--- Inspect the climate types in South America more closely
# Re-run the PCA with three axes on the (non-log-transformed) climate and colour
# every site by its position in that 3-D climate space.

include("prepare_data.jl")
bioclim_sa = prepare_environment(datadir)
pr2, load = do_pca(bioclim_sa, asm; naxes = 3)

biplot(pr2[:,1], pr2[:,2], load[:,1:2])

ct = pr2 .- minimum(pr2, dims = 1)
ct ./= maximum(ct)
ct .+= 0.5 .* (1 .- maximum(ct, dims = 1))

cols = [RGB(sl...) for sl in eachrow(ct)]
Plots.scatter(Tuple.(cellindices(asm)), color = cols, msw = 0, ms = 2, aspect_ratio = 1, yflip = true, size = (800, 1100), legend = false)

Plots.scatter(pca1, pca2, color = cols, msw = 0, ms = 1, aspect_ratio = 1, size = (600, 600), legend = false)
save_figures && savefig("figures/climate_colors2.png")
