#
#
#

using NearestNeighbors

include("binning.jl")
include("geometry.jl")

"""
    all_triplets(xyz, w, rmin, rmax, Nbin, dp, Np)

    # Input
    - xyz::Array{Float, 2}: x, y, z coordinates of the first point (3xN).
    - w::Array{Float, 1}: weights for each point.
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Flaot: Number of line-of-sight bins

    # Output
    - hist::Array{Float, 5}: the weighted histogram of triplets.
"""
function all_triplets(xyz, w, rmin, rmax, Nbin, dp, Np; logr, logp)
    hist = zeros(Nbin, Nbin, Nbin, Np, Np)
    xyz_tree = KDTree(xyz)
    # First galaxy
    for i1 in 1:size(xyz)[2]
        idxs2 = inrange(xyz_tree, xyz[:,i1], rmax)
        for i2 in idxs2
            # Ignore neighbours that are before the first point (don't double count)
            if i2 < i1 continue end
            # For each second neighbour look for all third neighbours
            idxs3 = inrange(xyz_tree, xyz[:,i2], rmax)
            for i3 in idxs3
                # Ignore neighbours that are before the second point (don't double count)
                # Ignore the ones that are also not in range of the first point
                if i3 < i2 || !(i3 in idxs2) continue end
                # histogram here
                p12, r12 = par_perp_distance(xyz[:,i1], xyz[:,i2])
                p23, r23 = par_perp_distance(xyz[:,i2], xyz[:,i3])
                _, r31 = par_perp_distance(xyz[:,i3], xyz[:,i1])
                indexes = histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np, logr=logr, logp=logp)
                if nothing in indexes continue end
                ir12, ir23, ir31, ip12, ip23 = indexes
                hist[ir12, ir23, ir31, ip12, ip23] += w[i1]*w[i2]*w[i3]
            end
        end
    end
    return hist
end