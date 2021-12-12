#
#
#
using DelimitedFiles

include("../src/triple_loop.jl")
include("../src/geometry.jl")

function radecz_to_hist(ra, dec, z, w, ofile)
    xyz = radecz_to_xyz(ra, dec, z)
    hist = all_triplets(xyz, w, 0.1, 20, 10, 1, 40, logr=true, logp=false)
    hist_red = fold_histogram(hist)
    writedlm(ofile, hist_red)
end
