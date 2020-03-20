using DrWatson
@quickactivate "QMC"

using QMC

# main.jl
#
# A projector QMC program for the TFIM

using Random
Random.seed!(1234)

using ProgressMeter

using Measurements
using Statistics
using FFTW

using DelimitedFiles
using JLD2
using Printf
using Lattices

using DataStructures
using ArgParse


s = ArgParseSettings()


@add_arg_table! s begin
    "dims"
        help = "The dimensions of the square lattice"
        required = true
        arg_type = Int
        nargs = '+'
    "--periodic", "-p"
        help = "Periodic BCs"
        action = :store_true
    "--field"
        help = "Strength of the transverse field"
        arg_type = Float64
        default = 1.0
    "--interaction", "-J"
        help = "Strength of the interaction"
        arg_type = Float64
        default = 1.0
    "-M"
        help = "Half-size of the operator list"
        arg_type = Int
        default = 1000
    "--measurements", "-n"
        help = "Number of samples to record"
        arg_type = Int
        default = 100_000
    "--skip", "-s"
        help = "Number of MC steps to perform between each measurement"
        arg_type = Int
        default = 0
end

parsed_args = parse_args(ARGS, s)

println(parsed_args)

PBC = parsed_args["periodic"]
h = parsed_args["field"]
J_ = parsed_args["interaction"]

Dim = length(parsed_args["dims"])
nX = parsed_args["dims"]

BC = PBC ? Periodic : Fixed
BC_name = PBC ? "PBC" : "OBC"

if Dim == 1
    nX = nX[1]
    lattice = Chain(nX; boundary = BC)
elseif Dim == 2
    lattice = Square(nX...; boundary = BC)
else
    error("Unsupported number of dimensions")
end

# MC parameters
M = parsed_args["M"] # length of the operator_list is 2M
MCS = parsed_args["measurements"] # the number of samples to record
EQ_MCS = div(MCS, 10)
skip = parsed_args["skip"]  # number of MC steps to perform between each msmt

root = "../data/$(Dim)D/$(nX)/$(BC_name)/J$(J_)/h$(h)/skip$(skip)/"
mkpath(root)
samples_file = "$(root)samples.txt"
info_file = "$(root)info.txt"
qmc_state_file = "$(root)state.jld2"

# *******  Globals

H = TFIM(lattice, h, J_)
qmc_state = BinaryQMCState(H, M)

## FINITE-BETA
# beta = 20.0
# @showprogress "Warm up..." for i in 1:EQ_MCS
#     mc_step_beta!(qmc_state, H, beta; eq = true)

#     println(length(qmc_state.operator_list))
# end


measurements = zeros(Int, MCS, nspins(H))
mags = zeros(MCS)
ns = zeros(MCS)

# @showprogress "MCMC...   " for i in 1:MCS # Monte Carlo Steps
#     ns[i] = mc_step_beta!(qmc_state, H, beta) do cluster_data, qmc_state, H
#         #mags[i] = magnetization(qmc_state.left_config)
#         #spin_prop = sample(qmc_state)
#         spin_prop = qmc_state.left_config
#         mags[i] = magnetization(spin_prop)
#     end
# end

#include("error.jl")
#mag_sqr = mean_and_stderr(x -> x^2, mags)
#println(mag_sqr)

@showprogress "Warm up..." for i in 1:EQ_MCS
    mc_step!(qmc_state, H)
end

@showprogress "MCMC...   " for i in 1:MCS # Monte Carlo Production Steps
   mc_step!(qmc_state, H) do cluster_data, qmc_state, H
       spin_prop = sample(qmc_state)
       measurements[i, :] = spin_prop

       ns[i] = num_single_site_diag(qmc_state.operator_list)
       mags[i] = magnetization(spin_prop)
   end

   for _ in 1:skip
       mc_step!(qmc_state, H)
   end
end

# open(samples_file, "w") do io
#     writedlm(io, measurements, " ")
# end

# @time @save qmc_state_file qmc_state

energy(H::TFIM) = n -> ((H.h != 0) ? (-H.h * ((1.0 / n) - 1)) : 0.0) + H.J * (nbonds(H) / nspins(H))
mag = mean_and_stderr(mags)
abs_mag = mean_and_stderr(abs, mags)
mag_sqr = mean_and_stderr(x -> x^2, mags)
# energy = -(mean_and_stderr(ns) / beta)/nspins(H) + H.J * (nbonds(H) / nspins(H)) + H.h
# heat_cap = jackknife(ns, ns .^ 2) do ns, ns_sqrd
#     (ns_sqrd - ns^2 - ns) / nspins(H)
# end

# heat_cap = (mean(abs2, ns) - mean(ns)^2 - mean(ns)) / nspins(H)


# println((mag, abs_mag, mag_sqr, energy, heat_cap))


@time energy_ = jackknife(energy(H), ns)
@time corr_time = correlation_time(mags .^ 2)

# println()

open(info_file, "w") do file_io
    streams = [Base.stdout, file_io]

    for io in streams
        @printf(io, "⟨M⟩   = % .16f +/- %.16f\n", mag.val, mag.err)
        @printf(io, "⟨|M|⟩ = % .16f +/- %.16f\n", abs_mag.val, abs_mag.err)
        @printf(io, "⟨M^2⟩ = % .16f +/- %.16f\n", mag_sqr.val, mag_sqr.err)
        @printf(io, "⟨E⟩   = % .16f +/- %.16f\n\n", energy_.val, energy_.err)

        println(io, "Correlation time: $(corr_time)\n")

        println(io, "Operator list length: $(2*M)")
        println(io, "Number of MC measurements: $(MCS)")
        println(io, "Number of equilibration steps: $(EQ_MCS)")
        println(io, "Number of skips between measurements: $(skip)\n")

        println(io, "Dim = $(Dim), nX = $(nX), PBC = $(PBC), h = $(h), J = $(J_)\n")

        println(io, "Samples outputted to file: $(samples_file)")
    end
end