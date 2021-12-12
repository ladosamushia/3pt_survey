#
#
#

using Random

include("../src/triple_loop.jl")

Ngal = 2000000
xyz = rand(3, Ngal)*2000
w = ones(Ngal)

hist = all_triplets(xyz, w, 0.1, 20, 10, 1, 40, logr=true, logp=false)
sum(hist)