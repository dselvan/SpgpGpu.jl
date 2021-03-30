using SpgpGpu
using Test

@testset "SpgpGpu.jl" begin
    x = 2
    y = 2
    @test SpgpGpu.sum_values(x, y) == 4
end
