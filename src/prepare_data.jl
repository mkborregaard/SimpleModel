# Build the model's data object: a raster-backed SpatialEcology `Assemblage` of
# South American birds, carrying the two environmental PCA axes as site
# statistics. 

using Rasters
using RasterDataSources
using ArchGDAL
using MultivariateStats
using DataFrames
using GeoInterface; const GI = GeoInterface
using Extents
using Shapefile
using StatsBase
using LinearAlgebra
using JLD2
using SpatialEcology

# Find the right rotation for the PCA - copied from factorloadingmatrices (varimax)
function vmax(A::AbstractMatrix{TA}; gamma = 1.0, minit = 20, maxit = 1000,
    reltol = 1e-12) where TA
    d, m = size(A)
    m == 1 && return A
    # Warm up step: start with a good initial orthogonal matrix T by SVD and QR
    T = Matrix{TA}(I, m, m)
    B = A * T
    L,_,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
    T = L * M'
    if norm(T - Matrix{TA}(I, m, m)) < reltol
        T = Matrix(qr(randn(m, m)).Q)
        B = A * T
    end

    # Iteration step: get better T to maximize the objective (as described in Factor Analysis book)
    D = 0
    for k in 1:maxit
        Dold = D
        L,s,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
        T = L * M'
        D = sum(s)
        B = A * T
        if (abs(D - Dold)/D < reltol) && k >= minit
            break
        end
    end
    return T
end

# Download the bioclim variables for south america
function prepare_environment(datadir)
    ENV["RASTERDATASOURCES_PATH"] = joinpath(datadir, "Rasterdatasources")
    bioclim = RasterStack(CHELSA{BioClim}; lazy=true, version = 1)
    bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
    Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 20)
end

# Fit a PCA model to the climate and extract the two primary components, one
# value per site of `asm` (in site order, via the assemblage's cell indices).
function do_pca(bioclim_sa, asm; naxes = 2)
    ci = cellindices(asm)
    bc = permutedims(bioclim_sa, (X, Y))                 # canonical order, as the sites
    big_mat = permutedims(reduce(hcat, maplayers(A -> zscore(A[ci]), bc)))
    model = fit(PCA, big_mat; maxoutdim = naxes)
    pred = MultivariateStats.transform(model, big_mat)
    vm = vmax(loadings(model))
    pred2 = pred' * vm # the minus here and below are just a transformation to have high values top right
    -pred2, -(loadings(model) * vm)
end

# Load the bird shapefiles and pick the ones in South America
function loadranges(data::String, batches::Int, mask, datadir)
    shapefiles = [joinpath(datadir, data, "batch_$i.shp") for i in 1:batches]
    reduce(vcat, map(shapefiles) do sf
        df = DataFrame(Shapefile.Table(sf))
        filter(df) do row
            ext = GI.calc_extent(GI.trait(row.geometry), row.geometry)
            Extents.intersects(ext, Extents.extent(mask)) && row.origin == 1 && row.seasonal in 1:2 && row.presence in 1:4
        end
    end)
end

# a function to rasterize a species by name
function get_speciesmask(name, geoms, mask)
    ret = reduce(.|, map(findall(==(name), geoms.sci_name)) do i
        boolmask(geoms.geometry[i]; to = mask, boundary = :touches) .& mask
    end)
    rebuild(ret; name = name)
end

function prepare_data(datadir; doplots = false)
    ## Get the environmental data
    bioclim_sa = prepare_environment(datadir)
    sa_mask = boolmask(bioclim_sa.bio15)
    doplots && Plots.plot(bioclim_sa.bio1)

    # log transform the precipitation layers
    for prec in 12:19
        lay = Symbol("bio$prec")
        bioclim_sa[lay][sa_mask] .= log.(bioclim_sa[lay][sa_mask] .+ 1)
    end

    ## Get the bird data (this takes time) and build the assemblage
    sa_geoms = loadranges("Birds", 5, sa_mask, datadir)
    allspecies = collect(skipmissing(unique(sa_geoms.sci_name)))
    allranges = RasterSeries([get_speciesmask(name, sa_geoms, sa_mask) for name in allspecies], (; name = allspecies))
    asm = Assemblage(allranges, sa_mask)

    ## Attach the two environmental PCA axes as site statistics
    pcamat, loads = do_pca(bioclim_sa, asm)
    addsitestats!(asm, pcamat[:, 1], :pca1)
    addsitestats!(asm, pcamat[:, 2], :pca2)
    doplots && biplot(asm[!, :pca1], asm[!, :pca2], loads, string.(names(bioclim_sa)))

    JLD2.save(joinpath(datadir, "processed_assemblage.jld2"), "asm", asm)
    asm
end
