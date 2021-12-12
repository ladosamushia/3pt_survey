#
#
#

using Random

include("../src/triple_loop.jl")
include("../src/geometry.jl")

Ngal = 2000000

ra = rand(Ngal)*π/2
dec = -π/3 .+ rand(Ngal)*2*π/3
z = 0.5 .+ rand(Ngal)
w = ones(Ngal)

xyz = radecz_to_xyz(ra, dec, z)
hist = all_triplets(xyz, w, 0.1, 20, 10, 1, 40, logr=true, logp=false)
sum(hist)