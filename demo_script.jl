using SpgpGpu
using Statistics

x = parse.(Float64, readlines("test_data/train_inputs"))
y = parse.(Float64, readlines("test_data/train_outputs"))
xt = parse.(Float64, readlines("test_data/test_inputs"))
me_y = mean(y);
y0 = y .- me_y

n_observ, n_dim = size(x, 1), size(x, 2)

hyp = [1.2106, -0.51968, -2.7338]
xb = [1.1462, 0.11905, 5.3042, 1.441, 3.1853,
    2.89, 1.4412, 1.4414, 1.1465, 5.6914,
    3.1857, 1.6291, 1.4429, 1.6265, 4.3954,
    1.4429, 1.4428, 1.4416, 3.3293, 1.6156]

# PREDICTION
μ0, σ² = predict(y0, x, xb, xt, hyp)
μ = μ0 .+ me_y
σ² = σ² .+ exp(hyp[end])
