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
