using CUDA
using KernelAbstractions, CUDAKernels

function dist(x0, x1)
    D = zeros(length(x0), length(x1))
    @tullio D[i,j] = x0[i] - x1[j];
    
    return D
end

function dist2(x0, x1)
    n0 = length(x0);
    n1 = length(x1);
    D = repeat(x0,1,n1)-repeat(x1',n0,1);
    
    return D
end