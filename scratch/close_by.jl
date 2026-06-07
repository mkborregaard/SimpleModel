
gridcell = [(50*Int(floor(i/10)) + Int(floor(j/10))) for i in axes(env.mask,1), j in axes(env.mask,2) if env.mask[i,j]]

sp = rand(spec.names)
ran = spec.ranges[At(sp)][env.mask]

function plot_gc(gc::Int)
    @show gc
    a = gridcell .== gc
    plot_var(a, env)
end

function plot_var(var, env)
    Plots.scatter(env.pca1, env.pca2, c = :black,  msw = 0, ms = 0.8, legend = false)
    Plots.scatter!(env.pca1[var], env.pca2[var], c = :red,  msw = 0, ms = 2)
end

plot_gc(rand(gridcell))


# one interesting example
gc = 1415

# which species has this grid cell?
a = gridcell .== gc
sp = rand(spec.ranges)
while count(sp[env.mask] .& a) == 0
    sp = rand(spec.ranges)
end

