using Test

include("../src/triple_loop.jl")
include("../src/binning.jl")

Ngal = 100
L = 10
xyz = rand(3, Ngal)*L
w = rand(1, Ngal)*L
rmin = 1.0
rmax = 10
Nbin = 10
dp = 1
Np = 10
hist = all_triplets(xyz, w, rmin, rmax, Nbin, dp, Np; logr=true, logp=false)
hist_direct = all_triplets_direct(xyz, w, rmin, rmax, Nbin, dp, Np; logr=true, logp=false)
# comparison of non-folded histograms may fail because of the order of galaxies
# in the triple does not have to be the same
hist_f = fold_histogram(hist)
hist_f_direct = fold_histogram(hist_direct)
@test isapprox(hist_f, hist_f_direct)
