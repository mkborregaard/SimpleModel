using Rasters, Plots, ArchGDAL, RasterDataSources, Statistics

datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
ENV["RASTERDATASOURCES_PATH"] = joinpath(datadir, "Rasterdatasources")
bioclim = RasterStack(CHELSA{BioClim}; lazy=true, version = 1)
bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
bioclim_sa = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 4)
bs2 = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 20)

sa_mask = boolmask(bioclim_sa.bio15)

t1 = bioclim_sa.bio1[boolmask(bioclim_sa.bio15)]
p1 = bioclim_sa.bio13[boolmask(bioclim_sa.bio15)]
t2 = bs2.bio1[boolmask(bioclim_sa.bio15)]
p2 = bs2.bio13[boolmask(bioclim_sa.bio15)]


GLMakie.scatter(t1[sa_mask], p1[sa_mask])

lookup(t1, X)


sortperm(t1[sa_mask])


getcd(fineval, coarsevec) = searchsortedlast(>(fineval), coarsevec)

# consider what to do with missing values

