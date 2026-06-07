using Random
using RecipesBase
using Statistics
using LinearAlgebra
using StatsBase: weights, fit, Histogram
import GeometryBasics

# A basic Ellipse struct
struct Ellipse
    center_x::Float64
    center_y::Float64
    length::Float64
    width::Float64
    angle::Float64 # Given in radians
end

function distance(point::Tuple, el::Ellipse)
    cosa = cos(el.angle)
    sina = sin(el.angle)
    rel_x = first(point) - el.center_x
    rel_y = last(point) - el.center_y
    a = (cosa * rel_x + sina * rel_y)^2 / el.length^2
    b = (sina * rel_x + cosa * rel_y)^2 / el.width^2
    a+b
end

in_ellipse(point, el::Ellipse) = distance(point, el) <= 1

rescale(x, lo, hi) = clamp(muladd(x, hi - lo, lo), lo, hi)

import GeometryBasics
function GeometryBasics.coordinates(el::Ellipse, nvertices = 100)
    xs = range(-el.length, el.length, length = nvertices ÷ 2)
    ys = [sqrt(1 - (x/el.length)^2) * el.width for x in xs]
    xs = [xs; reverse(xs)[2:end]; first(xs)]
    ys = [ys; reverse(-ys)[2:end]; first(ys)]
    x_ret = xs .* cos.(el.angle) - ys .* sin.(el.angle)
    y_ret = xs .* sin.(el.angle) + ys .* cos.(el.angle)
    Point2.(zip(x_ret .+ el.center_x, y_ret .+ el.center_y))
end

RecipesBase.@recipe f(el::Ellipse; nvertices = 100) = GeometryBasics.coordinates(el, nvertices)

truncnormrand() = clamp(randn(), -2, 2) / 6 + 0.5
truncrand() = clamp(rand(), 0.3, 0.99)

import Random.rand

# a random ellipse around a given center
function Random.rand(::Type{Ellipse}, x::Number, y::Number; 
    area = 1, lengthfun = truncrand) 
    a = lengthfun() * sqrt(pi)
    b = 1 / (pi * a)
    a, b = extrema((a,b))
    Ellipse(x, y, a * sqrt(area), b * sqrt(area), rand()π)
end

# pick the center randomly between some limits
Random.rand(::Type{Ellipse}; xlims=(0,1), ylims=(0,1), area = 1, lengthfun = truncrand) = 
    rand(Ellipse, rescale(rand(), xlims...), rescale(rand(), ylims...); area, lengthfun)

GeometryBasics.area(el::Ellipse) = el.length * el.width * π

# possibly use a covariance matrix weighted
# by the 1 / number of point occurrences in
# the same pca grid cell?

function StatsBase.fit(::Type{Ellipse}, xs, ys, sigma = 2; weight = ones(length(xs)))
    evals, evecs = eigen(cov([xs ys], weights(weight)))
    a, b = sigma .* sqrt.(evals)
    angle = atan(evecs[2,1] / evecs[1,1])
    Ellipse(mean(xs, weights(weight)), mean(ys, weights(weight)), a, b, angle)
end



