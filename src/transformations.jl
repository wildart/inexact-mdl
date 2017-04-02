module DataTransform

    import StatsBase: fit, mean_and_std

    export ZScoreTransform, UnitRangeTransform, fit, transform, transform!, reconstruct, reconstruct!

    ### Transformtion
    abstract AbstractDataTransform

    # apply the transform
    transform!{T<:AbstractFloat, S<:AbstractDataTransform}(t::S, x::DenseArray{T,1}) = transform!(x, t, x)
    transform!{T<:AbstractFloat, S<:AbstractDataTransform}(t::S, x::DenseArray{T,2}) = transform!(x, t, x)

    transform{T<:Real, S<:AbstractDataTransform}(t::S, x::DenseArray{T,1}) = transform!(Array(Float64, size(x)), t, x)
    transform{T<:Real, S<:AbstractDataTransform}(t::S, x::DenseArray{T,2}) = transform!(Array(Float64, size(x)), t, x)

    # reconstruct the original data from transformed values
    reconstruct!{T<:AbstractFloat, S<:AbstractDataTransform}(t::S, x::DenseArray{T,1}) = reconstruct!(x, t, x)
    reconstruct!{T<:AbstractFloat, S<:AbstractDataTransform}(t::S, x::DenseArray{T,2}) = reconstruct!(x, t, x)

    reconstruct{T<:Real, S<:AbstractDataTransform}(t::S, y::DenseArray{T,1}) = reconstruct!(Array(Float64, size(y)), t, y)
    reconstruct{T<:Real, S<:AbstractDataTransform}(t::S, y::DenseArray{T,2}) = reconstruct!(Array(Float64, size(y)), t, y)

    # Z-score transformation
    immutable ZScoreTransform <: AbstractDataTransform
        dim::Int
        mean::Vector{Float64}
        scale::Vector{Float64}

        ZScoreTransform() = new(0, Float64[], Float64[])

        function ZScoreTransform(d::Int, m::Vector{Float64}, s::Vector{Float64})
            lenm = length(m)
            lens = length(s)
            lenm == d || lenm == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
            lens == d || lens == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
            new(d, m, s)
        end
    end

    indim(t::ZScoreTransform) = t.dim
    outdim(t::ZScoreTransform) = t.dim

    # fit a z-score transform
    function fit{T<:Real}(::Type{ZScoreTransform}, X::DenseArray{T,2}; center::Bool=true, scale::Bool=true)
        d, n = size(X)
        n >= 2 || error("X must contain at least two columns.")

        m, s = mean_and_std(X, 2)

        return ZScoreTransform(d, (center ? convert(Vector{Float64}, vec(m)) : Array(Float64, 0)),
                                  (scale ? convert(Vector{Float64}, vec(s)) : Array(Float64, 0)))
    end


    function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,1}, t::ZScoreTransform, x::DenseArray{XT,1})
        d = t.dim
        length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if !is(x, y)
                    copy!(y, x)
                end
            else
                for i = 1:d
                    @inbounds y[i] = x[i] / s[i]
                end
            end
        else
            if isempty(s)
                for i = 1:d
                    @inbounds y[i] = x[i] - m[i]
                end
            else
                for i = 1:d
                    @inbounds y[i] = (x[i] - m[i]) / s[i]
                end
            end
        end
        return y
    end

    function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,2}, t::ZScoreTransform, x::DenseArray{XT,2})
        d = t.dim
        size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(x,2)
        size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if !is(x, y)
                    copy!(y, x)
                end
            else
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    # xj = x[:, j]
                    # yj = y[:, j]
                    for i = 1:d
                        @inbounds yj[i] = xj[i] / s[i]
                    end
                end
            end
        else
            if isempty(s)
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    # xj = x[:, j]
                    # yj = y[:, j]
                    for i = 1:d
                        @inbounds yj[i] = xj[i] - m[i]
                    end
                end
            else
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    # xj = x[:, j]
                    # yj = y[:, j]
                    for i = 1:d
                        @inbounds yj[i] = (xj[i] - m[i]) / s[i]
                    end
                end
            end
        end
        return y
    end

    function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,1}, t::ZScoreTransform, y::DenseArray{YT,1})
        d = t.dim
        length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if !is(y, x)
                    copy!(x, y)
                end
            else
                for i = 1:d
                    @inbounds x[i] = y[i] * s[i]
                end
            end
        else
            if isempty(s)
                for i = 1:d
                    @inbounds x[i] = y[i] + m[i]
                end
            else
                for i = 1:d
                    @inbounds x[i] = y[i] * s[i] + m[i]
                end
            end
        end
        return x
    end

    function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,2}, t::ZScoreTransform, y::DenseArray{YT,2})
        d = t.dim
        size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(y,2)
        size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        m = t.mean
        s = t.scale

        if isempty(m)
            if isempty(s)
                if !is(y, x)
                    copy!(x, y)
                end
            else
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    for i = 1:d
                        @inbounds xj[i] = yj[i] * s[i]
                    end
                end
            end
        else
            if isempty(s)
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    for i = 1:d
                        @inbounds xj[i] = yj[i] + m[i]
                    end
                end
            else
                for j = 1:n
                    xj = view(x, :, j)
                    yj = view(y, :, j)
                    for i = 1:d
                        @inbounds xj[i] = yj[i] * s[i] + m[i]
                    end
                end
            end
        end
        return x
    end

    # UnitRangeTransform normalization

    immutable UnitRangeTransform  <: AbstractDataTransform
        dim::Int
        unit::Bool
        min::Vector{Float64}
        scale::Vector{Float64}

        UnitRangeTransform() = new(0, false, Float64[], Float64[])

        function UnitRangeTransform(d::Int, unit::Bool, min::Vector{Float64}, max::Vector{Float64})
            lenmin = length(min)
            lenmax = length(max)
            lenmin == d || lenmin == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
            lenmax == d || lenmax == 0 || throw(DimensionMismatch("Inconsistent dimensions."))
            new(d, unit, min, max)
        end
    end

    indim(t::UnitRangeTransform) = t.dim
    outdim(t::UnitRangeTransform) = t.dim

    function fit{T<:Real}(::Type{UnitRangeTransform}, X::DenseArray{T,2}; unit::Bool=true)
        d, n = size(X)

        tmin = Array(Float64, d)
        tmax = Array(Float64, d)
        copy!(tmin, X[:, 1])
        copy!(tmax, X[:, 1])
        for j = 2:n
            @inbounds for i = 1:d
                if X[i, j] < tmin[i]
                    tmin[i] = X[i, j]
                elseif X[i, j] > tmax[i]
                    tmax[i] = X[i, j]
                end
            end
        end
        for i = 1:d
            @inbounds tmax[i] = 1.0 / (tmax[i] - tmin[i])
        end
        return UnitRangeTransform(d, unit, tmin, tmax)
    end

    function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,1}, t::UnitRangeTransform, x::DenseArray{XT,1})
        d = t.dim
        length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            for i = 1:d
                @inbounds y[i] = (x[i] - tmin[i]) * tscale[i]
            end
        else
            for i = 1:d
                @inbounds y[i] = x[i] * tscale[i]
            end
        end
        return y
    end

    function transform!{YT<:Real,XT<:Real}(y::DenseArray{YT,2}, t::UnitRangeTransform, x::DenseArray{XT,2})
        d = t.dim
        size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(x,2)
        size(y,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds yj[i, j] = (xj[i, j] - tmin[i]) * tscale[i]
                    # @inbounds y[i, j] = (x[i, j] - tmin[i]) * tscale[i]
                end
            end
        else
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds yj[i] = xj[i, j] * tscale[i]
                end
            end
        end
        return y
    end

    function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,1}, t::UnitRangeTransform, y::DenseArray{YT,1})
        d = t.dim
        length(x) == length(y) == d || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            for i = 1:d
                @inbounds x[i] = y[i] / tscale[i] +  tmin[i]
            end
        else
            for i = 1:d
                @inbounds x[i] = y[i] / tscale[i]
            end
        end
        return x
    end

    function reconstruct!{YT<:Real,XT<:Real}(x::DenseArray{XT,2}, t::UnitRangeTransform, y::DenseArray{YT,2})
        d = t.dim
        size(x,1) == size(y,1) == d || throw(DimensionMismatch("Inconsistent dimensions."))
        n = size(y,2)
        size(x,2) == n || throw(DimensionMismatch("Inconsistent dimensions."))

        tmin = t.min
        tscale = t.scale

        if t.unit
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds xj[i] = yj[i] / tscale[i] + tmin[i]
                end
            end
        else
            for j = 1:n
                xj = view(x, :, j)
                yj = view(y, :, j)
                for i = 1:d
                    @inbounds xj[i] = yj[i] / tscale[i]
                end
            end
        end
        return x
    end

end


# using DataTransform
#=
using Base.Test

@testset "Data transformations" begin

    X = rand(5, 8)

    @testset "Z-score transformation" begin
        t = fit(ZScoreTransform, X; center=false, scale=false)
        Y = transform(t, X)
        @test isa(t, DataTransform.AbstractDataTransform)
        @test isempty(t.mean)
        @test isempty(t.scale)
        @test isequal(X, Y)
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X

        t = fit(ZScoreTransform, X; center=false)
        Y = transform(t, X)
        @test isa(t, DataTransform.AbstractDataTransform)
        @test isempty(t.mean)
        @test length(t.scale) == 5
        @test_approx_eq Y X ./ std(X, 2)
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X

        t = fit(ZScoreTransform, X; scale=false)
        Y = transform(t, X)
        @test isa(t, DataTransform.AbstractDataTransform)
        @test length(t.mean) == 5
        @test isempty(t.scale)
        @test_approx_eq Y X .- mean(X, 2)
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X

        t = fit(ZScoreTransform, X)
        Y = transform(t, X)
        @test isa(t, DataTransform.AbstractDataTransform)
        @test length(t.mean) == 5
        @test length(t.scale) == 5
        @test_approx_eq Y (X .- mean(X, 2)) ./ std(X, 2)
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X
    end

    @testset "Unit range transformation" begin
        t = fit(UnitRangeTransform, X)
        Y = transform(t, X)
        @test length(t.min) == 5
        @test length(t.scale) == 5
        @test_approx_eq Y (X .- minimum(X, 2)) ./ (maximum(X, 2) .- minimum(X, 2))
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X

        t = fit(UnitRangeTransform, X; unit=false)
        Y = transform(t, X)
        @test length(t.min) == 5
        @test length(t.scale) == 5
        @test_approx_eq Y X ./ (maximum(X, 2) .- minimum(X, 2))
        @test_approx_eq transform(t, X[:,1]) Y[:,1]
        @test_approx_eq reconstruct(t, Y[:,1]) X[:,1]
        @test_approx_eq reconstruct(t, Y) X
    end
end
=#