#
#
#

using Random

include("../src/triple_loop.jl")

xyz = rand(3, 1000000)*2000
w = ones(100)

hist = all_triplets(xyz, w, 0.1, 20, 20, 1, 40, logr=true, logp=false)