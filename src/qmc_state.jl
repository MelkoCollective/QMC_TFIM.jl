function init_op_list(length)
    operator_list::Vector{NTuple{2,Int}} = [(0, 0) for _ in 1:length]
    return operator_list
end

function resize_op_list!(operator_list::Vector{NTuple{2, Int}}, new_size::Int)
    len = length(operator_list)

    if len < new_size
        tail = init_op_list(new_size - len)
        append!(operator_list, tail)
    end
end

abstract type AbstractQMCState{D,N} end

abstract type AbstractGroundState{D,N} <: AbstractQMCState{D,N} end
abstract type AbstractThermalState{D,N} <: AbstractQMCState{D,N} end

struct BinaryGroundState{N} <: AbstractGroundState{2,N}
    left_config::BitArray{N}
    right_config::BitArray{N}
    propagated_config::BitArray{N}

    operator_list::Vector{NTuple{2,Int}}

    linked_list::Vector{Int}
    leg_types::BitVector
    associates::Vector{NTuple{3,Int}}

    first::Vector{Int}
end


function BinaryGroundState(H::Hamiltonian{2,N}, M::Int) where N
    operator_list = init_op_list(2*M)

    len = 2*nspins(H) + 4*length(operator_list)
    linked_list = zeros(Int, len)
    leg_types = falses(len)
    associates = [(0, 0, 0) for _ in 1:len]

    first = zeros(Int, nspins(H))

    BinaryGroundState{N}(zero(H), zero(H), zero(H),
                         operator_list,
                         linked_list, leg_types, associates,
                         first)
end


function BinaryGroundState(left_config::BitArray{N}, right_config::BitArray{N}, operator_list::Vector{NTuple{2,Int}}) where N
    @assert left_config !== right_config "left_config and right_config can't be the same array!"

    len = 2*length(left_config) + 4*length(operator_list)
    linked_list = zeros(Int, len)
    leg_types = falses(len)
    associates = [(0, 0, 0) for _ in 1:len]

    first = zeros(Int, length(left_config))

    BinaryGroundState{N}(left_config, right_config, copy(left_config),
                         operator_list,
                         linked_list, leg_types, associates,
                         first)
end


struct BinaryThermalState{N} <: AbstractThermalState{2,N}
    left_config::BitArray{N}
    right_config::BitArray{N}
    propagated_config::BitArray{N}

    operator_list::Vector{NTuple{2,Int}}

    linked_list::Vector{Int}
    leg_types::BitVector
    associates::Vector{NTuple{3,Int}}

    first::Vector{Int}
    last::Vector{Int}
end


function BinaryThermalState(H::Hamiltonian{2,N}, cutoff::Int) where N
    operator_list = init_op_list(cutoff)

    len = 4*length(operator_list)
    linked_list = zeros(Int, len)
    leg_types = falses(len)
    associates = [(0, 0, 0) for _ in 1:len]

    first = zeros(Int, nspins(H))
    last = zeros(Int, nspins(H))

    BinaryThermalState{N}(zero(H), zero(H), zero(H),
                          operator_list,
                          linked_list, leg_types, associates,
                          first, last)
end


function BinaryThermalState(left_config::BitArray{N}, right_config::BitArray{N}, operator_list::Vector{NTuple{2,Int}}) where N
    @assert left_config !== right_config "left_config and right_config can't be the same array!"

    len = 4*length(operator_list)
    linked_list = zeros(Int, len)
    leg_types = falses(len)
    associates = [(0, 0, 0) for _ in 1:len]

    first = zeros(Int, length(left_config))
    last = copy(first)

    BinaryThermalState{N}(left_config, right_config, copy(left_config),
                          operator_list,
                          linked_list, leg_types, associates,
                          first, last)
end


const BinaryQMCState{N} = Union{BinaryGroundState{N}, BinaryThermalState{N}}


struct ClusterData
    # linked_list::Vector{Int}
    # leg_types::BitVector
    # associates::Vector{NTuple{3,Int}}
    # first::Vector{Int}
    last::Union{Vector{Int}, Nothing}
end