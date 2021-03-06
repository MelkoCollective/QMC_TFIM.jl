using BenchmarkTools
using QMC

SUITE = BenchmarkGroup()

SUITE["probability_vector"] = BenchmarkGroup()

for n in 10:10:100
    p = rand(1:10, n)

    p_vec = ProbabilityVector(p)
    p_heap = ProbabilityHeap(p)

    SUITE["probability_vector"][n] = BenchmarkGroup()
    SUITE["probability_vector"][n]["vector"] = @benchmarkable rand($p_vec)
    SUITE["probability_vector"][n]["heap"] = @benchmarkable rand($p_heap)
end


SUITE["TFIM_groundstate"] = BenchmarkGroup()
SUITE["TFIM_thermalstate"] = BenchmarkGroup()

for M = 200:200:1000
    bonds, Ns, Nb = lattice_bond_spins(10)
    H = TFIM(bonds, 1, Ns, Nb, 1.0, 1.0)

    ###########################################################################

    SUITE["TFIM_groundstate"][M] = BenchmarkGroup()
    groundstate = BinaryGroundState(H, M)

    SUITE["TFIM_groundstate"][M]["diagonal_update"] =
        @benchmarkable QMC.diagonal_update!($groundstate, $H)
    SUITE["TFIM_groundstate"][M]["linked_list_update"] =
        @benchmarkable QMC.linked_list_update($groundstate, $H)
    SUITE["TFIM_groundstate"][M]["cluster_update"] =
        @benchmarkable(QMC.cluster_update!(cd, $groundstate, $H),
                       setup=(cd = QMC.linked_list_update($groundstate, $H)))

    SUITE["TFIM_groundstate"][M]["mc_step"] = @benchmarkable QMC.mc_step!($groundstate, $H)

    ###########################################################################

    SUITE["TFIM_thermalstate"][M] = BenchmarkGroup()
    beta = 10.0
    thermalstate = BinaryThermalState(H, M)

    SUITE["TFIM_thermalstate"][M]["diagonal_update"] =
        @benchmarkable QMC.diagonal_update_beta!($thermalstate, $H, $beta)
    SUITE["TFIM_thermalstate"][M]["linked_list_update"] =
        @benchmarkable QMC.linked_list_update_beta($thermalstate, $H)
    SUITE["TFIM_thermalstate"][M]["cluster_update"] =
        @benchmarkable(QMC.cluster_update_beta!(cd, $thermalstate, $H),
                       setup=(cd = QMC.linked_list_update_beta($thermalstate, $H)))

    SUITE["TFIM_thermalstate"][M]["mc_step"] = @benchmarkable QMC.mc_step_beta!($thermalstate, $H, $beta)
end
