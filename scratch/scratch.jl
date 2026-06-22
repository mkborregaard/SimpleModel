
# Test the ellipse code #TODO move to plotting and script
using Plots
pts = vec(collect(Iterators.product(1:5:1000, 1:5:1000)))
el = rand(Ellipse; xlims=(1, 1000), ylims=(1, 1000), lengthlims=(1,400), widthlims=(1,400))
Plots.scatter(pts; aspect_ratio = 1, marker_z = [in_ellipse(x, el) for x in pts], msw = 0, ms = 1)

el = rand(Ellipse; xlims = (-5,10), ylims = (-5,5), lengthlims = 1.5, widthlims = 1.5)
els = [in_ellipse(pt, el) for pt in zip(pcas.pca1, pcas.pca2)]
h = Makie.scatter(collect(zip(pca1, pca2)); markersize = 0.1, color = :grey)
collect(zip(pcas.pca1, pcas.pca2))[els]
Makie.scatter!(collect(zip(pcas.pca1, pcas.pca2))[els]; markersize = 0.1, color = :red)

# Project niche pca elipses back to real space as masks
Rasters.rplot(els)

as, bs, es = Float64[], Float64[], Float64[]
for i in 1:100
    e = rand(Ellipse)
    a, b = e.length, e.width
    push!(es, a/b)
    push!(as, a)
    push!(bs, b)
end
Plots.scatter(as, bs)
Plots.histogram(es)



# plot a species on the density of points in climate space
background = fit(Histogram, (pca1, pca2),  (-10:0.3:10, -10:0.3:10))
spec = rand(allspecies)
x,y = get_climate(allranges[At(spec)], env)
sp = fit(Histogram, (x,y),  (-10:0.3:10, -10:0.3:10))
Plots.plot(Plots.plot(sp, aspect_ratio = 1), Plots.plot(background, aspect_ratio = 1), layout = 2)

# relative occupancy of climate space
rel = sp.weights ./ background.weights
replace!(rel, 0 => NaN)
Plots.heatmap(rel')




# not really necessary
#using GeometryBasics
#simplified = GeometryOps.simplify(sa_geoms[1:100]; tol=0.0000001);
#GeoInterface.convert(GeometryBasics, simplified)
#rast = rasterize(count, sa_geoms; to=sa_mask, threaded=false)
#Plots.plot(rast)

