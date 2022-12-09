using SpgpGpu
using LinearAlgebra, Statistics, StatsBase
using CUDA
using BenchmarkTools

use_gpu = false; # didn't quite get this working in time :/
train_first = false;

x = parse.(Float64, readlines("test_data/train_inputs"))
y = parse.(Float64, readlines("test_data/train_outputs"))
xtest = parse.(Float64, readlines("test_data/test_inputs"))
me_y = mean(y);
y0 = y .- me_y

N, dim = size(x, 1), size(x, 2)

if train_first
    M = 20
    idx = sortperm(vec(rand(N, 1)))
    idx = idx[1:M]
    xb_init = x[idx, :]

    hyp_init = zeros(dim + 2, 1)
    hyp_init[1:dim, 1] .= -2 * log((maximum(x) - minimum(x))' / 2)
    hyp_init[dim+1, 1] = log(var(y0, 1))
    hyp_init[dim+2, 1] = log(var(y0, 1) / 4)

    w_init = [reshape(xb_init, M * dim, 1); hyp_init]
    w, f = minimize(w_init, "spgp_lik", -200, y0, x, M)
    xb = reshape(w[1:M*dim, 1], M, dim)
    hyp = w[M*dim+1:end, 1]
else

    hyp = [1.2106, -0.51968, -2.7338]
    xb = [1.1462, 0.11905, 5.3042, 1.441, 3.1853,
        2.89, 1.4412, 1.4414, 1.1465, 5.6914,
        3.1857, 1.6291, 1.4429, 1.6265, 4.3954,
        1.4429, 1.4428, 1.4416, 3.3293, 1.6156]
end

if use_gpu
    y0 = cu(y0)
    x = cu(x)
    xb = cu(xb)
    xtest = cu(xtest)
    hyp = cu(hyp)
    del = cu(1e-6)
end

# PREDICTION
μ0, σ² = spgp_pred(y0, x, xb, xtest, hyp);
μ = μ0 .+ me_y;
σ² = σ² .+ exp(hyp[end]);
mu=μ;
@btime spgp_pred($y0, $x, $xb, $xtest, $hyp)
