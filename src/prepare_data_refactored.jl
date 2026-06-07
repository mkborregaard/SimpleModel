# Replaces the bottom half of prepare_data.jl (from "Get the bird data" onward).
# The environment prep and PCA functions are unchanged.
#
# Key change: instead of returning Dict(:spec => Species(...), :env => Environment(...))
# we return Dict(:asm => Assemblage(...), :env => EnvSummary(...))

function prepare_data(datadir; doplots = false)
    ## Environment (unchanged) -------------------------------------------------
    bioclim_sa = prepare_environment(datadir)
    sa_mask    = boolmask(bioclim_sa.bio15)

    doplots && Plots.plot(bioclim_sa.bio1)

    for prec in 12:19
        lay = Symbol("bio$prec")
        bioclim_sa[lay][sa_mask] .= log.(bioclim_sa[lay][sa_mask] .+ 1)
    end

    pcamat, loads = do_pca(bioclim_sa, sa_mask)
    pca1, pca2 = pcamat[:,1], pcamat[:,2]
    doplots && biplot(pca1, pca2, loads, string.(names(bioclim_sa)))

    pca_maps = RasterStack((pca1 = do_map(pca1, sa_mask), pca2 = do_map(pca2, sa_mask)))
    doplots && Plots.plot(pca_maps)

    ## Species (builds Assemblage) ---------------------------------------------
    sa_geoms   = loadranges("Birds", 5, sa_mask, datadir)
    allspecies = collect(skipmissing(unique(sa_geoms.sci_name)))
    allranges  = RasterSeries(
        [get_speciesmask(name, sa_geoms, sa_mask) for name in allspecies],
        (; name = allspecies)
    )

    ## Per-site metadata as a DataFrame ----------------------------------------
    inds = collect(Iterators.product(1:size(sa_mask,1), 1:size(sa_mask,2)))[sa_mask]
    ss   = DataFrame(pca1 = pca1, pca2 = pca2, inds = inds)

    ## Build the Assemblage (new) -----------------------------------------------
    asm = Assemblage(allranges, sa_mask; sitestats = ss)

    ## Domain-level summaries (small, not per-site) ----------------------------
    bbox  = (extrema(pca1), extrema(pca2))
    cc    = concave_hull(ConcaveHull.KDTree(vcat(pca1', pca2'), reorder = false), 40)
    chull = GI.Polygon([GI.LinearRing([cc.vertices; [first(cc.vertices)]])])
    env   = EnvSummary(bbox, chull)

    ## Save and return ----------------------------------------------------------
    obj = Dict("asm" => asm, "env" => env, "pca_maps" => pca_maps)
    JLD2.save(joinpath(datadir, "processed_objects.jld2"), obj)
    obj
end
