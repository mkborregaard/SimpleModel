using Rasters
using GLMakie
using Stencils
using GeometryOps

el = rand(ellipses)
els = [in_ellipse(pt, el) for pt in zip(pca1, pca2)]
elmap = do_map(els, sa_mask; missingval = false)

# Polygonize the masks
# This needs the improve_polygonize branch of GeometryOps
centers = map(lookup(els)) do lookup
    parent(Rasters.shiftlocus(Rasters.Center(), lookup))
end
polygons = GeometryOps.polygonize(centers..., els)
Rasters.rplot(els)
GLMakie.poly!(polygons)

for p in polygons
    GLMakie.poly!(p)
end

# Check the polygon areas - potential range sizes
res = map(abs ∘ step, centers)
min_viable_area = 20 # Maybe this is really small? # MKB: Actually a lot of the ranges are really really small
too_small_removed = filter(polygons) do polygon
    count(boolmask(polygon; res)) > min_viable_area
end 
Rasters.rplot(els)
Makie.poly!(too_small_removed)

# Join areas withing some allowed dispersal radius
dilate(a) = Stencils.mapstencil(Stencils.Circle{4}(), a) do hood
    any(hood)
end
dilated_polygons = GeometryOps.polygonize(centers..., dilate(els))
Rasters.rplot(els)
Makie.poly!(dilated_polygons)

multi_polygons = GI.MultiPolygon.(map(dilated_polygons) do dp
    filter(polygons) do p
        LibGEOS.within(p, dp)
    end 
end |> filter(x -> length(x) > 0));

# Plot habitat groups that may be the same species
GLMakie.poly(dilated_polygons)
for mp in multi_polygons
    Makie.poly!(mp)
end

p = Rasters.rplot(pcas.pca1; transparency=true)
GLMakie.poly!.(p.axis, multi_polygons)
GLMakie.save("species.png", p)
p = GLMakie.scatter( collect(zip(pca1, pca2)); markersize = 0.1, color = :grey)
collect(zip(pcas.pca1, pcas.pca2))[els]
GLMakie.scatter!(collect(zip(pcas.pca1, pcas.pca2))[els]; markersize = 0.1, color = :red)
GLMakie.save("pca.png", p)

