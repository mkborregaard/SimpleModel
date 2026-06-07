using Rasters
using GeoInterface

struct Environment{N <: Number, S <: RasterStack, R <: Raster{Bool}, P <: GeoInterface.Wrappers.Polygon}
    pca1::Vector{N}                         # the first environmental pca axis of two (rotated to both hold equal variance)
    pca2::Vector{N}                         # the second environmental pca axis
    pca_maps::S                             # the pca values as Raster objects for easy plotting
    mask::R                                 # a Raster{Bool} outlining the domain (e.g. South America)
    inds::Vector{Tuple{Int, Int}}           # a Vector that translates indices (into e.g. the pca vectors) to indices into the raster maps
    bbox::Tuple{Tuple{N, N}, Tuple{N, N}}   # bounding box for the points in pca space
    chull::P                                # concave hull for the points in pca space
end

struct Species{R <: RasterSeries}
    ranges::R                               # a Series of Bool Rasters, each one specifying the distribution of one species at the resolution of `mask`
    names::Vector{String}                   # names to index into the ranges
end

function Base.show(io::IO, e::Environment)
    println(io, "Environment with $(length(e.pca1)) sites on a $(size(e.mask, 1))x$(size(e.mask, 2)) domain")
end