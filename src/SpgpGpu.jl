module SpgpGpu

using LinearAlgebra, Tullio, Printf

export spgp_pred
include("spgp_pred.jl")

export kernel
include("kernel.jl")

export dist, dist2
include("dist.jl")

end
