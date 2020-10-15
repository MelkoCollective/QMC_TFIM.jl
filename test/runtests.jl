using Test, QMC
using Statistics

@testset "Small Vector of size $n" for n in 1:10
    p = rand(1:10, n)
    mean_ = p' * collect(1:size(p,1)) / sum(p)
    var_ = (p' * (collect(1:size(p,1)) .^ 2) / sum(p)) - mean_^2

    @testset "ProbabilityVector" begin
        pvec = ProbabilityVector(p)
        X = [rand(pvec) for _ in 1:1_000_000]

        @test isapprox(mean_, mean(X), atol=0.05)
        @test isapprox(var_, var(X), atol=0.05)

        # X = [QMC.rand_bisect(pvec) for _ in 1:1_000_000]
        # @test isapprox(mean_, mean(X), atol=0.05)
        # @test isapprox(var_, var(X), atol=0.05)
    end

    @testset "ProbabilityHeap" begin
        pvec = ProbabilityHeap(p)
        X = [rand(pvec) for _ in 1:1_000_000]

        @test isapprox(mean_, mean(X), atol=0.05)
        @test isapprox(var_, var(X), atol=0.05)
    end
end
