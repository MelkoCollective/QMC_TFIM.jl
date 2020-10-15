# (-3, i) = long field
# (-2, i) = transverse field off-diag
# (-1, i) = transverse field diag
# (0, 0) = id
# (i, j) = diag bond op

function make_prob_vector(J::AbstractMatrix{T}, hx::AbstractVector{T}) where T
    ops = Vector{NTuple{2, Int}}(undef, 0)
    p = Vector{T}(undef, 0)

    k = 0
    for i in eachindex(hx)
        if hx[i] != 0
            push!(ops, (-1, i))
            push!(p, h)
        end
    end

    # only take J_ij terms from upper diagonal since it's symmetric
    for j in axes(J, 2), i in axes(J, 1)
        if i < j
            if J[i, j] != 0
                push!(ops, op)
                push!(p, 2*abs(J))
            end
        end
    end

    return ops, p
end

function make_prob_vector(bond_spins::Vector{NTuple{2,Int}}, Ns::Int, J::T, h::T) where T
    ops = Vector{NTuple{2, Int}}(undef, 0)
    p = Vector{T}(undef, 0)

    if !iszero(h)
        for i in 1:Ns
            push!(ops, (-1, i))
            push!(p, h)
        end
    end

    if !iszero(J)
        for op in bond_spins
            push!(ops, op)
            push!(p, 2*abs(J))
        end
    end

    return ops, p
end


###################


abstract type AbstractOperatorSampler{N, T, P <: AbstractProbabilityVector{T}} end
firstindex(::AbstractOperatorSampler) = 1
lastindex(os::AbstractOperatorSampler) = length(os)


struct OperatorSampler{N, T, P} <: AbstractOperatorSampler{N, T, P}
    operators::Vector{NTuple{N, Int}}
    pvec::P
end


# samples operators using a heap
# based on a blog post by Tim Vieira
# https://timvieira.github.io/blog/post/2016/11/21/heaps-for-incremental-computation/
function OperatorSampler(operators::Vector{NTuple{N, Int}}, p::Vector{T}) where {N, T <: Real}
    @assert length(operators) == length(p) "Given vectors must have the same length!"
    pvec = probability_vector(p)
    return OperatorSampler{N, T, typeof(pvec)}(operators, pvec)
end
rand(os::OperatorSampler{N}) where N = @inbounds os.operators[rand(os.pvec)]
@inline length(os::OperatorSampler) = length(os.operators)




function cluster_probs_vec(operators::Vector{NTuple{N, Int}}, p::AbstractVector{T}) where {N, T <: Real}
    perm = sortperm(p)
    p = p[perm]
    operators = operators[perm]

    uniq_p = T[]
    uniq_ops = Vector{Vector{NTuple{N, Int}}}(undef, 0)

    for i in axes(p, 1)
        if length(uniq_p) == 0
            push!(uniq_p, p[i])
            push!(uniq_ops, [operators[i]])
            continue
        end

        if !(last(uniq_p) ≈ p[i])
            push!(uniq_p, p[i])
            push!(uniq_ops, [operators[i]])
        else
            push!(uniq_ops[end], operators[i])
        end
    end

    # rescale uniq_p
    for i in eachindex(uniq_p, uniq_ops)
        uniq_p[i] *= length(uniq_ops[i])
    end
    return uniq_ops, uniq_p
end



struct HierarchicalOperatorSampler{N, T, P} <: AbstractOperatorSampler{N, T, P}
    operator_bins::Vector{Vector{NTuple{N, Int}}}
    pvec::P
end

function HierarchicalOperatorSampler(operators::Vector{NTuple{N, Int}}, p::AbstractVector{T}) where {N, T <: Real}
    @assert length(operators) == length(p) "Given vectors must have the same length!"
    operator_bins, p = cluster_probs_vec(operators, p)
    pvec = probability_vector(p)
    return HierarchicalOperatorSampler{N, T, typeof(pvec)}(operator_bins, pvec)
end
length(os::HierarchicalOperatorSampler) = sum(length, os.operator_bins)

function rand(os::HierarchicalOperatorSampler{N})::NTuple{N, Int} where N
    @inbounds ops_list = os.operator_bins[rand(os.pvec)]
    l = rand(1:length(ops_list))

    return @inbounds ops_list[l]
end
