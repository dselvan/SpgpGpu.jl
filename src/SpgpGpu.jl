module SpgpGpu

using LinearAlgebra

export predict
include("predict.jl")

export kernel
include("kernel.jl")

end
