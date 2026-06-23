# Turn a fitted niche ellipse into a contiguous geographic range, by flood-filling
# outward from a random occupied cell across the cells whose climate falls inside
# the ellipse.

using Rasters
using SpatialEcology: to_raster
import GeometryBasics

const nbh = Tuple((x, y) for x in -1:1, y in -1:1 if !(x == y == 0))

on_domain(pt, domain) = min(pt...) > 0 && first(pt) <= size(domain, 1) && last(pt) <= size(domain, 2) && domain[pt...]

function add_point!(georange, potentials, pt, domain)
    georange[pt...] = true
    for nb in nbh
        newpt = nb .+ pt
        if on_domain(newpt, domain) && !georange[newpt...]
            push!(potentials, newpt)
        end
    end
end

function growrange(start, domain)
    on_domain(start, domain) || error("start point not on domain")
    georange = fill(false, dims(domain), missingval = false)
    potentials = Set([start])
    failsafe = 0; max_iter = prod(size(domain))
    while !isempty(potentials)
        (failsafe += 1) > max_iter && error("stuck")
        pt = pop!(potentials)
        add_point!(georange, potentials, pt, domain)
    end
    georange
end

function random_point_on_ellipse(el::Ellipse, x, y; maxiter = 1e6)
    iter = 0
    while (iter += 1) < maxiter
        pt = rand(1:length(x))
        in_ellipse((x[pt], y[pt]), el) && return pt
    end
    error("Did not find a point on the ellipse in $maxiter tries")
end
random_point_on_ellipse(el::Ellipse, asm; maxiter = 1e6) =
    random_point_on_ellipse(el, asm[!, :pca1], asm[!, :pca2]; maxiter)

# A Bool raster marking every site whose climate lies inside the ellipse.
map_ellipse(el::Ellipse, asm) =
    to_raster([in_ellipse(pt, el) for pt in zip(asm[!, :pca1], asm[!, :pca2])], asm)

function make_continuous_range(el, asm)
    domain = map_ellipse(el, asm)
    GeometryBasics.area(el) == 0 && return domain
    i = random_point_on_ellipse(el, asm[!, :pca1], asm[!, :pca2])
    pt = Tuple(cellindices(asm)[i])
    growrange(pt, domain)
end
