abstract type AbstractProbabilityVector{T <: Real} <: AbstractVector{T} end
Base.@propagate_inbounds getindex(pvec::AbstractProbabilityVector, i) = getweight(pvec, i)
Base.@propagate_inbounds setindex!(pvec::AbstractProbabilityVector, w, i) = setweight!(pvec, w, i)
Base.@propagate_inbounds setprobability!(pvec::AbstractProbabilityVector{T}, p::T, i::Int) where T = setweight!(pvec, p*normalization(pvec), i)
Base.@propagate_inbounds getprobability(pvec::AbstractProbabilityVector, i) = getweight(pvec, i) / normalization(pvec)

firstindex(::AbstractProbabilityVector) = 1
lastindex(pvec::AbstractProbabilityVector) = length(pvec)
size(pvec::AbstractProbabilityVector) = (length(pvec),)

###############################################################################

# more efficient for small probability vectors
struct ProbabilityVector{T} <: AbstractProbabilityVector{T}
    p::Vector{T}
    cdf::Vector{T}

    function ProbabilityVector{T}(p::Vector{T}) where {T <: Real}
        if length(p) == 0
            throw(ArgumentError("probability vector must have non-zero length!"))
        end
        if any(x -> x < zero(T), p)
            throw(ArgumentError("weights must be non-negative!"))
        end
        cdf = cumsum(p)
        new{T}(p, cdf)
    end
end
ProbabilityVector(p::Vector{T}) where T = ProbabilityVector{T}(p)
@inline length(pvec::ProbabilityVector) = length(pvec.p)
normalization(pvec::ProbabilityVector) = @inbounds pvec.cdf[end]

function show(io::IO, pvec::ProbabilityVector{T}) where T
    r = repr(pvec.p; context=IOContext(io, :limit=>true))
    println(io, "ProbabilityVector{$T}($r)")
end


function rand(pvec::ProbabilityVector{T})::Int where T
    cdf = pvec.cdf
    r = rand() * cdf[end]

    for i in eachindex(cdf)
        @inbounds if r < cdf[i]
            return i
        end
    end
    return lastindex(cdf)
end

# TODO: fix bisection method
# function rand_bisect(pvec::ProbabilityVector{T})::Int where T
#     if length(pvec) < 5
#         return rand_linear(pvec)
#     else
#         cdf = pvec.cdf
#         # @show length(cdf)
#         r = rand()*cdf[end]
#         l, u = firstindex(pvec)-1, lastindex(pvec)-1

#         while l <= u
#             m = fld((l + u), 2)
#             cdfm = cdf[m+1]
#             if cdfm < r
#                 l = m + 1
#             elseif cdfm > r
#                 u = m - 1
#             else
#                 return m + 1
#             end
#         end
#         return l
#     end
# end

@inline function setweight!(pvec::ProbabilityVector{T}, w::T, i::Int) where T
    @boundscheck checkbounds(pvec, i)

    p, cdf = pvec.p, pvec.cdf
    @inbounds diff = w - p[i]
    @inbounds p[i] = w

    @simd for j in i:length(p)
        cdf[j] += diff
    end
end


@inline function getweight(pvec::ProbabilityVector, i::Int)
    @boundscheck checkbounds(pvec, i)
    return @inbounds pvec.p[i]
end

@inline function getweight(pvec::ProbabilityVector, r::AbstractArray{Int})
    @boundscheck checkbounds(pvec, r)
    return @inbounds pvec.p[r]
end


###############################################################################

# more efficient for larger probability vectors
struct ProbabilityHeap{T} <: AbstractProbabilityVector{T}
    prob_heap::Vector{T}
    length::Int

    function ProbabilityHeap{T}(p::Vector{T}) where {T <: Real}
        if length(p) == 0
            throw(ArgumentError("probability vector must have non-zero length!"))
        end
        if any(x -> x < zero(T), p)
            throw(ArgumentError("weights must be non-negative!"))
        end

        L = length(p)
        d = 2 ^ ceil(Int, log2(L))

        heap = zeros(T, 2*d)
        @inbounds @views heap[d : d + L - 1] = p

        @inbounds for i in (d-1):-1:1
            heap[i] = heap[2*i] + heap[2*i + 1]
        end

        new{T}(heap, L)
    end
end
ProbabilityHeap(p::Vector{T}) where T = ProbabilityHeap{T}(p)
@inline length(pvec::ProbabilityHeap) = pvec.length
normalization(pvec::ProbabilityHeap) = @inbounds pvec.prob_heap[1]

function show(io::IO, p::ProbabilityHeap{T}) where T
    heap, L = p.prob_heap, p.length

    d = 2 ^ ceil(Int, log2(L))

    pvec = heap[d : d + L - 1]
    r = repr(pvec; context=IOContext(io, :limit=>true))

    println(io, "ProbabilityHeap{$T}($r)")
end


function rand(pvec::ProbabilityHeap{T})::Int where T
    heap = pvec.prob_heap
    l = length(heap) ÷ 2
    @inbounds r = rand() * heap[1]

    i = 1
    while i < l
        i *= 2
        @inbounds left = heap[i]
        if r > left
            r -= left
            i += 1
        end
    end
    return (i - l + 1)
end

@inline function setweight!(pvec::ProbabilityHeap{T}, w::T, i::Int) where T
    @boundscheck checkbounds(pvec, i)

    heap = pvec.prob_heap
    j = (length(heap) ÷ 2) - 1 + i
    @inbounds heap[j] = w
    @inbounds while j > 1
        j ÷= 2
        heap[j] = heap[2*j] + heap[2*j + 1]
    end
end

@inline function getweight(pvec::ProbabilityHeap, i::Int)
    @boundscheck checkbounds(pvec, i)

    heap = pvec.prob_heap
    d = (length(heap) ÷ 2) - 1
    @inbounds x = heap[d + i]
    return x
end

@inline function getweight(pvec::ProbabilityHeap, r::AbstractArray{Int})
    @boundscheck checkbounds(pvec, r)

    heap = pvec.prob_heap
    d = (length(heap) ÷ 2) - 1
    @inbounds x = heap[d .+ r]
    return x
end

###############################################################################


const CUTOFF = 50

function probability_vector(p::Vector{T})::AbstractProbabilityVector{T} where T
    if length(p) < CUTOFF
        return ProbabilityVector(p)
    else
        return ProbabilityHeap(p)
    end
end
