"""
    predict(y, x, xb, xt, hyp, ∇)

    predictive distribution for SPGP given hyperparameters and pseudo-inputs

    y -- training targets (N x 1)
    x -- training inputs (N x dim)
    xb -- pseudo inputs (n x dim)
    xt -- test inputs (Nt x dim)
    hyp -- hyperparameters (including noise)
           for Gaussian covariance: (dim+2 x 1)
           
          hyp(1:dim) = log( b )
          hyp(dim+1) = log( c )
          hyp(dim+2) = log( sig )

         where cov = c * exp[-0.5 * sum_d b_d * (x_d - x'_d)^2] 
                         + sig * delta(x,x')

    ∇ -- OPTIONAL jitter (default 1e-6)
    μ -- predictive mean
    σ² -- predictive variance of latent function
          (add noise if full variance is required)
"""
function spgp_pred(y, x, xb, xt, hyp, ∇=1e-6)
    n_observ, n_dim = size(x, 1), size(x, 2)
    n_pseudo = size(xb, 1)
    n_train = size(xt, 1)
    
    σ = exp(hyp[end])

    K = kernel(xb, xb, hyp) + ∇ * I
    L = cholesky(Hermitian(K)).L
    K = kernel(xb, x, hyp)
    V = L \ K

    @inline function kdiag(x, hyp)
        c = exp(hyp[end - 1])
        Kd = fill(c, size(x, 1), 1)
        
        return Kd
    end

    ɛ = 1 .+ (kdiag(x, hyp) - transpose(sum(V.^2, dims=1))) / σ

    V = V ./ repeat(transpose(sqrt.(ɛ)), n_pseudo, 1);
    y = y ./ sqrt.(ɛ);
    Lm = transpose(cholesky(σ * I + V * transpose(V)).U);
    bet = Lm \ (V * y);

    K = kernel(xb, xt, hyp);
    lst = L \ K;
    lmst = Lm \ lst;
    μ = transpose((transpose(bet) * lmst));

    σ² = kdiag(xt, hyp) .- transpose(sum(lst.^2, dims=1)) .+ σ .* transpose(sum(lmst.^2, dims=1));

    return μ, σ²
end