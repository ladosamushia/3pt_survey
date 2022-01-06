using Test

include("../src/binning.jl")

rmin = 0.1
rmax = 20
Nbin = 20
dp = 1
Np = 40
indexes = histogram(5, 6, 7, 2.1, 4.9, rmin, rmax, Nbin, dp, Np, logr=true, logp=false)
@test indexes == (15, 16, 17, 43, 45)
indexes = histogram(0.5, 5, 5, -2.1, 4.9, rmin, rmax, Nbin, dp, Np, logr=true, logp=false)
@test indexes == (7, 15, 15, 38, 45)
indexes = histogram(0.5, 5, 5, -2.1, -4.9, rmin, rmax, Nbin, dp, Np, logr=true, logp=false)
@test indexes == (7, 15, 15, 38, 36)
indexes = histogram(5, 6, 21, 2.1, 4.9, rmin, rmax, Nbin, dp, Np, logr=true, logp=false)
@test indexes == nothing

ind = bin(0.5, 0.1, 20, 20, logscale=true)
@test ind == 7
ind = bin(-1, 0.1, 20, 20, logscale=true)
@test ind == 9
ind = bin(25, 0.1, 20, 20, logscale=true)
@test isnothing(ind)
ind = bin(0.01, 0.1, 20, 20, logscale=true)
@test isnothing(ind)
ind = bin(3.1, 1, 10, 11, logscale=false)
@test ind == 3

hist = ones(10,10,10,5,5)
hist_red = fold_histogram(hist)
for i in 1:size(hist_red)[2]
    if hist_red[1,i] == hist_red[2,i] == hist_red[3,i]
        @test hist_red[6,i] == 1
    elseif hist_red[1,i] == hist_red[2,i] || hist_red[2,i] == hist_red[3,i] || hist_red[3,i] == hist_red[1,i]
        @test hist_red[6,i] == 3
    else
        @test hist_red[6,i] == 6
    end
end
