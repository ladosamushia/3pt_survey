#
# logp has to be false in this version.
#

using NearestNeighbors

include("binning.jl")
include("geometry.jl")

"""
    all_triplets(xyz, w, rmin, rmax, Nbin, dp, Np; logr, logp)

    # Input
    - xyz::Array{Float, 2}: x, y, z coordinates of the first point (3xN).
    - w::Array{Float, 1}: weights for each point.
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Flaot: Number of line-of-sight bins
    - logr::Bool: Is the r binning logarithmic?
    - logp::Bool: Is the pi binning logarithmic?
    # Output
    - hist::Array{Float, 5}: the weighted histogram of triplets.
    
    logr, logp get passed to histogram function.
    Only for self-triplets.
"""
function all_triplets(xyz, w, rmin, rmax, Nbin, dp, Np; logr, logp)
    hist = zeros(Nbin, Nbin, Nbin, 2*Np, 2*Np)
    xyz_tree = KDTree(xyz)
    Rmax = sqrt(rmax^2 + (dp*Np)^2) + 10 # a small buffer just in case
    # First galaxy
    for i1 in 1:size(xyz)[2]
        idxs2 = inrange(xyz_tree, xyz[:,i1], Rmax)
        for i2 in idxs2
            # Ignore neighbours that are before the first point (don't double count)
            if i2 < i1 continue end
            # For each second neighbour look for all third neighbours
            idxs3 = inrange(xyz_tree, xyz[:,i2], Rmax)
            for i3 in idxs3
                # Ignore neighbours that are before the second point (don't double count)
                # Ignore the ones that are also not in range of the first point
                if i3 < i2 || !(i3 in idxs2) continue end
                # histogram here
                p12, r12 = par_perp_distance(xyz[:,i1], xyz[:,i2])
                p23, r23 = par_perp_distance(xyz[:,i2], xyz[:,i3])
                p31, r31 = par_perp_distance(xyz[:,i3], xyz[:,i1])
                indexes = histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np, logr=logr, logp=logp)
		if indexes == nothing continue end
                ir12, ir23, ir31, ip12, ip23 = indexes
                hist[ir12, ir23, ir31, ip12, ip23] += w[i1]*w[i2]*w[i3]
            end
        end
    end
    return hist
end


"""
    all_triplets(xyz1, xyz2, w1, w2, rmin, rmax, Nbin, dp, Np, logr, logp)

    # Input
    - xyz1::Array{Float, 2}: x, y, z coordinates of the first point (3xN) from the first set.
    - xyz2::Array{Float, 2}: x, y, z coordinates of the second point (3xN) from the second set.
    - w1::Array{Float, 1}: weights for each point from the first set.
    - w2::Array{Float, 1}: weights for each point from the second set.
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Float: Number of line-of-sight bins
    - logr::Bool: Is the r binning logarithmic?
    - logp::Bool: Is the pi binning logarithmic?
    # Output
    - hist::Array{Float, 5}: the weighted histogram of triplets.

    logr, logp are passed to histogram.
    This is for mixed triplets where xyz1/w1 are used twice and xyz2/w2 is used
    ones (e.g. for DDR counts, D=1, R=2).
"""
function all_triplets(xyz1, xyz2, w1, w2, rmin, rmax, Nbin, dp, Np; logr, logp)
    hist = zeros(Nbin, Nbin, Nbin, 2*Np, 2*Np)
    xyz1_tree = KDTree(xyz1)
    xyz2_tree = KDTree(xyz1)
    Rmax = sqrt(rmax^2 + (dp*Np)^2) + 10 # small offset just in case
    # First galaxy
    for i1 in 1:size(xyz1)[2]
        idxs2 = inrange(xyz1_tree, xyz1[:,i1], Rmax)
        for i2 in idxs2
            # Ignore neighbours that are before the first point (don't double count)
            if i2 < i1 continue end
            # For each second neighbour look for all third neighbours from the second set
            idxs32 = inrange(xyz2_tree, xyz1[:,i2], Rmax)
            idxs31 = inrange(xyz2_tree, xyz1[:,i1], Rmax)
	    # Only points that are correct distance from the both
	    idxs3 = union(idxs32, idxs31)
            for i3 in idxs3
                # histogram here
                p12, r12 = par_perp_distance(xyz1[:,i1], xyz1[:,i2])
                p23, r23 = par_perp_distance(xyz1[:,i2], xyz2[:,i3])
                p31, r31 = par_perp_distance(xyz2[:,i3], xyz1[:,i1])
                indexes = histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np, logr=logr, logp=logp)
                if indexes == nothing continue end
                ir12, ir23, ir31, ip12, ip23 = indexes
                hist[ir12, ir23, ir31, ip12, ip23] += w1[i1]*w1[i2]*w2[i3]
            end
        end
    end
    return hist
end

"""
    all_triplets_direct(xyz, w, rmin, rmax, Nbin, dp, Np; logr, logp)

    # Input
    - xyz::Array{Float, 2}: x, y, z coordinates of the first point (3xN).
    - w::Array{Float, 1}: weights for each point.
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Flaot: Number of line-of-sight bins
    - logr::Bool: Is the r binning logarithmic?
    - logp::Bool: Is the pi binning logarithmic?
    # Output
    - hist::Array{Float, 5}: the weighted histogram of triplets.
    
    logr, logp get passed to histogram function.
    Only for self-triplets.
    The same as all_triplets but uses the slow but exact triple loop.
"""
function all_triplets_direct(xyz, w, rmin, rmax, Nbin, dp, Np; logr, logp)
    hist = zeros(Nbin, Nbin, Nbin, 2*Np, 2*Np)
    Ngal = size(xyz)[2]
    # First galaxy
    for i1 in 1:Ngal
	# Second galaxy but don't double count
        for i2 in i1:Ngal
	    # Third galaxy but don't double count
	    for i3 in i2:Ngal
                # histogram here
                p12, r12 = par_perp_distance(xyz[:,i1], xyz[:,i2])
                p23, r23 = par_perp_distance(xyz[:,i2], xyz[:,i3])
                p31, r31 = par_perp_distance(xyz[:,i3], xyz[:,i1])
                indexes = histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np, logr=logr, logp=logp)
                if indexes == nothing continue end
                ir12, ir23, ir31, ip12, ip23 = indexes
                hist[ir12, ir23, ir31, ip12, ip23] += w[i1]*w[i2]*w[i3]
            end
        end
    end
    return hist
end

"""
    all_triplets_direct(xyz1, xyz2, w1, w2, rmin, rmax, Nbin, dp, Np, logr, logp)

    # Input
    - xyz1::Array{Float, 2}: x, y, z coordinates of the first point (3xN) from the first set.
    - xyz2::Array{Float, 2}: x, y, z coordinates of the second point (3xN) from the second set.
    - w1::Array{Float, 1}: weights for each point from the first set.
    - w2::Array{Float, 1}: weights for each point from the second set.
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Float: Number of line-of-sight bins
    - logr::Bool: Is the r binning logarithmic?
    - logp::Bool: Is the pi binning logarithmic?
    # Output
    - hist::Array{Float, 5}: the weighted histogram of triplets.

    logr, logp are passed to histogram.
    This is for mixed triplets where xyz1/w1 are used twice and xyz2/w2 is used
    ones (e.g. for DDR counts, D=1, R=2).
    Uses the slow explicit triple loop. Only for testing purposes.
"""
function all_triplets_direct(xyz1, xyz2, w1, w2, rmin, rmax, Nbin, dp, Nd; logr, logp)
    hist = zeros(Nbin, Nbin, Nbin, 2*Np, 2*Np)
    Ngal1 = size(xyz1)[2]
    Ngal2 = size(xyz2)[2]
    # First galaxy
    for i1 in 1:Ngal1
	# Second galaxy from the first set (don't double count)
        for i2 in i1:Ngal1
	    # All third galaxies from the second set
            for i3 in 1:Ngal2
                # histogram here
                p12, r12 = par_perp_distance(xyz1[:,i1], xyz1[:,i2])
                p23, r23 = par_perp_distance(xyz1[:,i2], xyz2[:,i3])
                p31, r31 = par_perp_distance(xyz2[:,i3], xyz1[:,i1])
                indexes = histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np, logr=logr, logp=logp)
                if indexes == nothing continue end
                ir12, ir23, ir31, ip12, ip23 = indexes
                hist[ir12, ir23, ir31, ip12, ip23] += w1[i1]*w1[i2]*w2[i3]
            end
        end
    end
    return hist
end

