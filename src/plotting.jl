using Plots
using GLMakie
using StatsBase
using GeoInterfaceMakie
using GeoInterfaceRecipes


# visualize a PCA biplot
const climvars = Dict(
    :bio1  => "mean annual air temperature",
    :bio2  => "mean diurnal air temperature range",
    :bio3  => "isothermality",
    :bio4  => "temperature seasonality",
    :bio5  => "mean daily maximum air temperature of the warmest month",
    :bio6  => "mean daily minimum air temperature of the coldest month",
    :bio7  => "annual range of air temperature",
    :bio8  => "mean daily mean air temperatures of the wettest quarter",
    :bio9  => "mean daily mean air temperatures of the driest quarter",
    :bio10  => "mean daily mean air temperatures of the warmest quarter",
    :bio11  => "mean daily mean air temperatures of the coldest quarter",
    :bio12  => "annual precipitation amount",
    :bio13  => "precipitation amount of the wettest month",
    :bio14  => "precipitation amount of the driest month",
    :bio15  => "precipitation seasonality",
    :bio16  => "mean monthly precipitation amount of the wettest quarter",
    :bio17  => "mean monthly precipitation amount of the driest quarter",
    :bio18  => "mean monthly precipitation amount of the warmest quarter",
    :bio19  => "mean monthly precipitation amount of the coldest quarter"
)

function biplot(pca1, pca2, loadings, labels = string.(keys(climvars)))
    p = GLMakie.scatter(collect(zip(pca1, pca2)); color=(:grey, 0.5))
    for i in 1:19
        GLMakie.lines!([(0,0), 5 .* (loadings[i,:]...,)], color = :red)
        GLMakie.text!((5 .* loadings[i,:])..., text = labels[i])
    end
    p
end

# This plots the distribution of a random species on the map
function plot_distribution(speciesname, spec::Species, env::Environment)
    Plots.plot(env.mask, fill = :grey)
    Plots.plot!(spec.ranges[At(speciesname)])
end

# Creates two side-by-side plots, one in geographic space, the other in climate space
# Shows the occurrence points of all the species defined by speciesnames
function plot_species(speciesnames, spec::Species, env::Environment)
    f = Figure(resolution = (1500, 700))
    a = Axis(f[1,1], aspect = DataAspect())
    b = Axis(f[1,2], aspect = DataAspect())
    Makie.plot!(a, env.mask, colormap = "Greys")
    Makie.scatter!(b, collect(zip(env.pca1, env.pca2)); markersize = 1, color=(:grey, 0.5))
    colsize!(f.layout, 1, Auto(0.73))
    rich = fill(NaN, dims(env.mask))

    for (i, species) in enumerate(speciesnames)
        species_mask = spec.ranges[At(species)]
        rich[species_mask] .= i
        Makie.scatter!(b, env.pca_maps.pca1[species_mask], env.pca_maps.pca2[species_mask]; 
            markersize = 2, label = species,
            colormap=(:tab10, 0.5), color=i, colorrange=(1, 10)
        )
    end

    Makie.plot!(a, rich, colormap=(:tab10, 0.5), colorrange=(1, 10))
  #  leg = Makie.Legend(f[1,3], a, framevisible = false, colorrange = (1,10)) # TODO: there was a problem with the legend code
  #  for i in 1:length(speciesnames)
  #      leg.entrygroups[][1][2][i].elements[1].markersize = 20
  #  end
    f
end


# Overplot some points on the pca background
function overplot_pca_space(points, color, env::Environment)
    p = Figure()
    a = Axis(p[1,1], aspect = DataAspect())
    Makie.scatter!(a, collect(zip(env.pca1, env.pca2)); markersize = 0.3, color=(:grey, 0.5))
    Makie.scatter!(a, points; markersize = 4, color = color, label = "rangesizes", colormap = Reverse(:Spectral))
    Makie.Colorbar(p[1,2], colormap = Reverse(:Spectral))
    p
end

# Overplot some points on the map of south america
function overplot_geo_space(points, col, env::Environment)
    f = Figure()
    a = Axis(f[1,1], aspect = DataAspect())
    Makie.plot!(a, env.mask, colormap = :Greys)
    Makie.scatter!(points, color = col)
    f
end

function plot_species_pca(speciesname, spec::Species, env::Environment, sigma = 2; weightmap = nothing)
    xs, ys = get_climate(speciesname, spec, env)
    length(xs) < 2 && return(Plots.scatter(env.pca1, env.pca2, ms = 1, mc = :grey, msw = 0, label = "", title = speciesname, aspect_ratio = 1))
    el = fit(Ellipse, xs, ys, sigma; weight = !isnothing(weightmap) ? weightmap[spec.ranges[At(speciesname)]] : ones(length(xs)))
    p = Plots.scatter(env.pca1, env.pca2, ms = 1, mc = :grey, msw = 0, label = "", title = speciesname, aspect_ratio = 1)
    Plots.scatter!(p, xs, ys, ms = 2, mc = :red, msw = 0, label = "occurrences")
    Plots.plot!(p, el, color = :blue, lw = 1, label = "niche")
    p
end

function plot_ellipse_patches(myel::Int, spec::Species, env::Environment, elsize = 2)
    Plots.plot(
        Plots.heatmap(first(env.pca_maps), title = "temp"),
        Plots.heatmap(last(env.pca_maps), title = "precip"),
        Plots.plot(spec.ranges[myel], title = "actual range"),
        plot_species_pca(spec.names[myel], spec, env, elsize; weightmap), 
        Plots.plot(map_ellipse(emp_ellipses[myel], env), title = "full ellipse range"),
        Plots.plot(make_continuous_range(emp_ellipses[myel], env), title = "sampled patch"),
        size = (800, 1000), layout = (3,2)
    )
end

