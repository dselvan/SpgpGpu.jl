function kernel(x1, x2, hyp)
    n1, dim = size(x1, 1), size(x1, 2)
    n2 = size(x2, 1)
    
    b = exp.(hyp[1:end - 2]) 
    c = exp.(hyp[end - 1]);

    x1 = x1 .* repeat(transpose(sqrt.(b)), n1, 1)
    x2 = x2 .* repeat(transpose(sqrt.(b)), n2, 1)
    
    K = (-2 * x1 * transpose(x2) +
        repeat(transpose(sum(x2 .* x2, dims=2)), n1, 1) +
        repeat(sum(x1 .* x1, dims=2), 1, n2)) 
    K = c * exp.(-0.5 * K)

    return K
end