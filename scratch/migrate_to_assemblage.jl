# One-time migration: convert the OLD cached object (`processed_objects.jld2`,
# a Dict of an `Environment` + a `Species`) into the new
# `processed_assemblage.jld2` (a single SpatialEcology `Assemblage`), so the
# refactored `script.jl` can load instantly instead of re-running `prepare_data`.
#
# Run once with:  julia --project=. scratch/migrate_to_assemblage.jl
# Safe to delete afterwards — it is only needed to reuse an existing old cache.

using Rasters
using JLD2
using DataFrames
using SpatialEcology

# Minimal stand-ins for the retired structs, so JLD2 can deserialize the old file.
struct Environment{N, S, R, P}
    pca1::Vector{N}; pca2::Vector{N}; pca_maps::S; mask::R
    inds::Vector{Tuple{Int, Int}}; bbox::Tuple{Tuple{N, N}, Tuple{N, N}}; chull::P
end
struct Species{R}
    ranges::R; names::Vector{String}
end

datadir = "/Users/cvg147/Dropbox/Arbejde/Data"

obj  = JLD2.load(joinpath(datadir, "processed_objects.jld2"))
env  = obj["env"]
spec = obj["spec"]

# The old per-site ordering (env.inds / env.pca1) matches the canonical (X, Y)
# cell order the Assemblage constructor produces, so the two PCA axes drop
# straight in as site statistics.
asm = Assemblage(spec.ranges, env.mask;
                 sitestats = DataFrame(pca1 = env.pca1, pca2 = env.pca2))

JLD2.save(joinpath(datadir, "processed_assemblage.jld2"), "asm", asm)
println("wrote processed_assemblage.jld2 — ",
        nspecies(asm), " species over ", nsites(asm), " sites")
