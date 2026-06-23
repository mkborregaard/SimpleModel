using Random
using RecipesBase
using Statistics
using LinearAlgebra
using StatsBase: weights, Histogram
import StatsBase: fit          # imported for extension (the Ellipse method below)
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
    b = (-sina * rel_x + cosa * rel_y)^2 / el.width^2
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
    GeometryBasics.Point2.(zip(x_ret .+ el.center_x, y_ret .+ el.center_y))
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

# ---------------------------------------------------------------------------
# Pluggable point-trust / quality-control strategies.
#
# Fitting a niche has two separable steps:
#   1. decide how much to *trust* each occurrence (it may be a false presence), then
#   2. estimate the ellipse shape and size from the trusted points.
# Step 1 is the open research question, so it is a strategy selected by multiple dispatch.
# To add a new SDM quality-control technique, define
#     struct MyQC <: TrustMethod ... end
#     trust(::MyQC, xs, ys; weight) = ...      # Vector{Float64} of per-point trust in [0,1]
# where 1 = fully trusted true presence and 0 = discarded. Hard filters return a 0/1 mask;
# soft methods may return fractional weights, which down-weight points in the fit.
# ---------------------------------------------------------------------------
abstract type TrustMethod end

# weighted (sigma = 1) ellipse: centre / orientation / aspect from the weighted moments.
# Returns a degenerate (zero) ellipse when there are too few trusted points or the
# covariance is not finite, so callers can detect and skip it via a zero length/width.
function shapefit(xs, ys, w)
    count(>(0), w) < 3 && return Ellipse(0, 0, 0, 0, 0)
    ww = weights(w)
    C = cov([xs ys], ww)
    all(isfinite, C) || return Ellipse(0, 0, 0, 0, 0)
    evals, evecs = eigen(C)
    Ellipse(mean(xs, ww), mean(ys, ww),
            sqrt(max(evals[1], 0)), sqrt(max(evals[2], 0)), atan(evecs[2,1], evecs[1,1]))
end

# No quality control: every occurrence is trusted.
struct TrustAll <: TrustMethod end
trust(::TrustAll, xs, ys; weight = ones(length(xs))) = ones(length(xs))

# Robust MCD-style concentration: from all points, repeatedly keep the `keep` fraction
# closest (Mahalanobis) to the running core and refit until the kept set converges.
# Iterative; resistant to outliers leveraging the shape. Returns a 0/1 mask.
Base.@kwdef struct MCDTrim <: TrustMethod
    keep::Float64 = 0.9
    maxiter::Int = 50
end
function trust(m::MCDTrim, xs, ys; weight = ones(length(xs)))
    n = length(xs)
    (m.keep >= 1 || n < 4) && return ones(n)
    h = max(3, round(Int, m.keep * n))
    keepset = trues(n); prev = Int[]
    for _ in 1:m.maxiter
        base = shapefit(xs, ys, weight .* keepset)
        (base.length == 0 || base.width == 0) && break
        order = sort(partialsortperm([distance((xs[i], ys[i]), base) for i in 1:n], 1:h))
        order == prev && break
        prev = order
        fill!(keepset, false); keepset[order] .= true
    end
    Float64.(keepset)
end

# Single-pass environmental-outlier filter: fit the niche to all points, then discard any
# point beyond the `p`-quantile of the chi-square (2 df) envelope. Parametric and fast;
# the 2-df quantile is closed form (-2 log(1-p)), so no extra dependency.
Base.@kwdef struct ChisqFilter <: TrustMethod
    p::Float64 = 0.95
end
function trust(m::ChisqFilter, xs, ys; weight = ones(length(xs)))
    base = shapefit(xs, ys, weight)
    base.length == 0 && return ones(length(xs))   # too few points to filter; trust all
    thr = -2 * log(1 - m.p)
    Float64.([distance((xs[i], ys[i]), base) <= thr for i in eachindex(xs)])
end

# Density filter: keep the `keep` fraction of points in the most densely occupied cells of
# climate space (2-D histogram at resolution `binsize`); discard records from sparse
# climates. Non-parametric, independent of any assumed niche shape.
Base.@kwdef struct DensityTrim <: TrustMethod
    keep::Float64 = 0.9
    binsize::Float64 = 0.2
end
function trust(m::DensityTrim, xs, ys; weight = ones(length(xs)))
    n = length(xs)
    (m.keep >= 1 || n < 4) && return ones(n)
    bins(v) = floor(minimum(v) - m.binsize):m.binsize:ceil(maximum(v) + m.binsize)
    hgram = fit(Histogram, (xs, ys), (bins(xs), bins(ys)))
    dens = [hgram.weights[searchsortedlast(hgram.edges[1], xs[i]),
                          searchsortedlast(hgram.edges[2], ys[i])] for i in 1:n]
    Float64.(dens .>= quantile(dens, 1 - m.keep))
end

# Fit an ellipse: trust the points via `method`, estimate the shape from the trusted
# (weighted) points, then size the boundary. With `cover` a fraction the boundary encloses
# that trust-weighted quantile of the points' Mahalanobis radii (cover = 1.0 -> all trusted
# points inside); with `cover === nothing` it sits at `sigma` core standard deviations.
function fit(::Type{Ellipse}, xs, ys, sigma = 2; weight = ones(length(xs)),
        method::TrustMethod = MCDTrim(), cover = 1.0)
    t = trust(method, xs, ys; weight)
    shape = shapefit(xs, ys, weight .* t)
    s = (cover === nothing || shape.length <= 0) ? sigma :
        weighted_quantile([sqrt(distance((xs[i], ys[i]), shape)) for i in eachindex(xs)], t, cover)
    Ellipse(shape.center_x, shape.center_y, s * shape.length, s * shape.width, shape.angle)
end

# weight-aware quantile: smallest x whose cumulative `w` mass reaches fraction `p` of the
# total. For a 0/1 trust mask this is the ordinary quantile over the trusted points.
function weighted_quantile(x, w, p)
    idx = sortperm(x)
    cw = cumsum(@view w[idx])
    total = last(cw)
    total <= 0 && return maximum(x)
    i = searchsortedfirst(cw, p * total)
    x[idx[clamp(i, 1, length(x))]]
end



