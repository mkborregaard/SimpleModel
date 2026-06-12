using Plots
using GLMakie
using Rasters
using JLD2

include("objects.jl")
include("plotting.jl")
include("ellipse.jl")
include("simplemodel.jl")
include("spread.jl")

save_figures = false

###--- First we load all the data into two objects. This takes a while if first time, 
# so we us

datadir = "/Users/cvg147/Dropbox/Arbejde/Data"
#datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
#datadir = "C:\\Users\\cvg147\\Dropbox\\Arbejde\\Data"
obj = try
    JLD2.load(joinpath(datadir, "processed_objects.jld2"))
catch
    include("prepare_data.jl")
    obj = prepare_data(datadir, doplots = true)
end

spec = obj["spec"] # see the "objects.jl" file for an explanation of these two objects
env = obj["env"]

###--- Exploratory data analysis

# Count total diversity
diversity = reduce(+, spec.ranges)
Plots.heatmap(diversity, color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/empirical_richness.png")


# Plot the diversity in climate space
f = Figure()
a, s = GLMakie.scatter(f[1,1], collect(zip(env.pca1, env.pca2)); markersize = 1, color = diversity[env.mask], colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2],s)
display(f)


# plot some random species and look at their distribution
randspecs = rand(spec.names, 5)
plot_species(randspecs, spec, env)


## how are ranges shaped?
cors, xrange, yrange = find_range_shapes(spec, env)
histogram(cors)
Plots.scatter(xrange, yrange)


###--- Patterns of range and niche size in climate space

# get environmental centroids and range sizes for all species
allcentroids = vec(get_centroid.(spec.ranges, (env,)))
rangesizes = vec(count.(spec.ranges))

# This figure shows the geographic range size of ranges on the environmental centroids (in environment spaecs)
overplot_pca_space(allcentroids, rangesizes, env)

# get hull area for all species
allhulls = [hullarea(sp, spec, env) for sp in spec.names]

# show the mean size of the climatic range in environmental space
overplot_pca_space(allcentroids, allhulls, env)

# get the geographical centroids for all species
gc = geocentroids.(spec.names, (spec,))

# plot the location of all the centroids on the map
overplot_geo_space(gc, rangesizes, env)

# divide the range sizes into quantiles for plotting
rangequants = asquantile(rangesizes, 4)

# and plot them
f = Figure()
aspect = DataAspect()
axs = [Axis(f[1,1]; aspect), Axis(f[2,1]; aspect), Axis(f[1,2]; aspect), Axis(f[2,2]; aspect)]
for i in 1:4
    inds = findall(==(i), rangequants)
    Makie.plot!(axs[i], env.mask, colormap = :Greys)
    Makie.scatter!(axs[i], gc[inds], markersize = 2)
end
f

# Now fit an ellipse to a species in pca space
plot_species_pca(rand(spec.names), spec, env, 2)

# Plot 16 random species with occurrences in pca space and fitted ellipses
p = Plots.plot([
    plot_species_pca(rand(spec.names), spec, env, 2) for i in 1:16]...
, size = (1200, 1200))
save_figures && savefig(p, "figures/16 species in pca space.png")

# Now control for the density of points in pca space by applying a 0.1 grid
# First map the options, e.g. binsizes from 0.1 to 0.4
map_binsize(binsize) = Plots.heatmap(1 ./ do_map(makeweights(env.pca1, env.pca2, binsize), env.mask), color = cgrad(:Spectral, rev = true), title = "binsize = $binsize")
Plots.plot([map_binsize(bs) for bs in 0.1:0.1:0.4]..., size = (800, 800))
save_figures && savefig("figures/binsizes.png")

# we conclude that we need one at 0.2
const weightmap = do_map(makeweights(env.pca1, env.pca2, 0.2), env.mask)

p = Plots.plot([
    plot_species_pca(rand(spec.names), spec, env, 1.5; weightmap) for i in 1:16]...
, size = (1200, 1200))
save_figures && savefig(p, "figures/16 species controlling for point density.png")

# fit elliptical niches for all species
emp_ellipses = [fitellipse(name, spec, env, 1.5) for name in spec.names]

# show patterns of ellipse area
ares = GeometryBasics.area.(emp_ellipses)
histogram(ares)
save_figures && savefig("figures/histogram of empirical ellipse areas.png")
Plots.scatter((el -> (el.center_x, el.center_y)).(emp_ellipses), marker_z = ares, ms = 3)
save_figures && savefig("figures/PCA centroids of empirical ellipses with area as color")


# repeat the plot with the 
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
el_emp_point = [count(el -> in_ellipse(pt, el), emp_ellipses) for pt in zip(env.pca1, env.pca2)]
Plots.plot(
    Plots.scatter(env.pca1, env.pca2, marker_z = el_emp_point, title = "fitted ellipse overlap"), 
    Plots.scatter(env.pca1, env.pca2, marker_z = diversity[env.mask], title = "empirical richness") 
)
save_figures && savefig("figures/empirical_ellipse_and_empirical_pca_richness.png")

  
# Create random ellipses with the empirical areas
rand_ellipses = [sample_ellipse(harea; env, on_real_point = true) for harea in ares]

# Show 50 random ellipses
p = Plots.scatter(env.pca1, env.pca2, mc = :grey, ms = 1, msw = 0, aspect_ratio = 1, label = "")
for el in rand(rand_ellipses, 50)
    Plots.plot!(p, el, label = "")
end
p
save_figures && savefig(p, "figures/50 random ellipses.png")

# plot the modelled and empirical richness
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
elpoint = [count(el -> in_ellipse(pt, el), rand_ellipses) for pt in zip(env.pca1, env.pca2)]
Plots.plot(
    Plots.scatter(env.pca1, env.pca2, marker_z = elpoint, title = "ellipse overlap"), 
    Plots.scatter(env.pca1, env.pca2, marker_z = diversity[env.mask], title = "empirical richness") 
)

save_figures && savefig("figures/modelled_ellipse_and_empirical_pca_richness.png")

Plots.heatmap(do_map(el_emp_point, env.mask), color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness_on_empirical_ellipses.png")

Plots.heatmap(do_map(elpoint, env.mask), color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness_on_random_ellipses.png")


###--- Let's create some random ranges and look at the patterns
# first with the spread model

# Let's look at a random ellipse just to know what is going on
plot_ellipse_patches(rand(eachindex(emp_ellipses)), spec, env, 1)

# range patches based on the empirical ellipses
model_ranges = RasterSeries([make_continuous_range(el, env) for el in emp_ellipses], (; name = spec.names))
newdiv = reduce(+, model_ranges)
Plots.heatmap(newdiv, color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness from patches in empirical ellipses.png")

# range patches based on randomly placed ellipses (with the right size)
rand_model_ranges = RasterSeries([make_continuous_range(el, env) for el in rand_ellipses], (; name = spec.names))
rand_div = reduce(+, rand_model_ranges)
Plots.heatmap(rand_div, color = cgrad(:Spectral, rev = true))
save_figures && savefig("figures/richness from patches in random ellipses.png")


# Find and plot richnes at the 1 degree lat/long scale
model_coarse = reduce(+, Rasters.aggregate.(any, model_ranges, 6))
emp_coarse = reduce(+, Rasters.aggregate.(any, spec.ranges, 6))

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


function occup_in_ellipse(el::Ellipse, env::Environment, sp::Raster{Bool})
    ins, outs = 0, 0
    for i in eachindex(env.pca1)
        if in_ellipse((env.pca1[i], env.pca2[i]), el)
            if sp[env.inds[i]...]
                ins += 1
            else
                outs += 1
            end
        end
    end
    ins, outs
end

allins, allouts = Int[], Int[]
for ind in eachindex(emp_ellipses)
    ins, outs = occup_in_ellipse(emp_ellipses[ind], env, spec.ranges[ind])
    push!(allins, ins)
    push!(allouts, outs)
end

Plots.scatter(log10.(GeometryBasics.area.(emp_ellipses) .+ 0.1), allins ./ (allins .+ allouts))
Plots.scatter(log10.(rangesizes .+ 1), log10.(allins ./ (allins .+ allouts)))
# it appears that rangesize is largely determined by how much of yuor niche you occupy
# or - niches are consistently overestimated (possibly more likely)
# some of the small-ranged species really occur in lots of regions - why is that?


Plots.scatter(log10.(GeometryBasics.area.(emp_ellipses) .+ 0.1), log10.(rangesizes))
# larger ranges have larger ellipses but not really that strong - I believe small ellipses are exaggerated

wd = findall(>(100), allins ./ allouts)
rangesizes[wd] # the second, number 103, is almost completely occupied. Let me take a look
plot_ellipse_patches(103, spec, env)
plot_ellipse_patches(3606, spec, env) # The Amazonian species tend to fully occupy their ellipse



# Inspect the climate types in South America more closely
# start with big_mat from prepare_environment

include("prepare_data.jl")
bioclim_sa = prepare_environment(datadir)
sa_mask = boolmask(bioclim_sa.bio15)
pr2, load = do_pca(bioclim_sa,sa_mask; naxes = 3)

biplot(pr2[:,1], pr2[:,2], load[:,1:2])

ct = pr2 .- minimum(pr2, dims = 1)
ct ./= maximum(ct)
ct .+= 0.5 .* (1 .- maximum(ct, dims = 1))

cols = [RGB(sl...) for sl in eachrow(ct)]
Plots.scatter(env.inds, color = cols, msw = 0, ms = 2, aspect_ratio = 1, yflip = 
true, size = (800, 1100), legend = false)

Plots.scatter(env.pca1, env.pca2, color = cols, msw = 0, ms = 1, aspect_ratio = 1, size = (600, 600), legend = false)

savefig("figures/climate_colors2.png")




