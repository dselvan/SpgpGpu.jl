using SpgpGpu
using Tullio 
using BenchmarkTools

x = vec(collect(1.0:10000.0));

# pure matlab imp in julia 
@btime y=dist2($x, $x); 

#using tullio einsum for a more julia way to do this
@btime y=dist($x, $x);

# now time for some GPU magic 
using CUDA, KernelAbstractions, CUDAKernels 

# define our function so the array preallocation can happen outside the function
gpu_dist(D,x0,x1) = @tullio D[i,j] = x0[i] - x1[j];
CUDA.allowscalar(false); # safety first, don't want to lock my laptop up during the demo
xg = CuArray(x); # copy the same array to gpu memory 
yg = CUDA.zeros(length(xg), length(xg));
@btime CUDA.@sync gpu_dist($yg, $xg, $xg)
