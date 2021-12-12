#
#
#


"""
    histogram(r12, r23, r31, p12, p23, hist, rmin, rmax, Nbin, dp, Np)

    # Input
    - r12::Float: distance between points 1 and 2
    - r23::Float: distance between points 2 and 3
    - r31::Float: distance between points 3 and 1
    - p12::Float: line-of-sight distance between 1 and 2
    - p23::Float: line-of-sight distance between 2 and 3
    - rmin::Float: minimum distance for binning
    - rmax::Float: maximum distance for binning
    - Nbin::Int: number of bins
    - dp::Float: line-of-sight bin width
    - Np::Flaot: Number of line-of-sight bins
    - logr::Bool: Is the r_per scale log?
    - logp::Bool: Is the r_par scale log?

    # Output
    - i1, i2, i3, i4, i5::Int: bin indexes

    The first three are binned in equally spaced log-bings.
    The last two in equally spaced linear bins (ignoring the sign).
"""
function histogram(r12, r23, r31, p12, p23, rmin, rmax, Nbin, dp, Np; logr, logp)
    i1 = bin(r12, rmin, rmax, Nbin, logscale=logr)
    i2 = bin(r23, rmin, rmax, Nbin, logscale=logr)
    i3 = bin(r31, rmin, rmax, Nbin, logscale=logr)
    i4 = bin(p12, 0, dp*Np, Np, logscale=logp)
    i5 = bin(p23, 0, dp*Np, Np, logscale=logp)
    return i1, i2, i3, i4, i5
end

"""
    bin(r, rmin, rmax, Nbin, logscale)

    # Input
    - r::Float: the distance
    - rmin::Float: where to start the binning
    - rmax::Float: where to end the binning
    - Nbin::Int: number of bins
    - logscale::Bool: is the binning on log scale?

    # Output
    - ind::Int: bin index for log binning

    Return nothing if out of the range.
"""
function bin(r, rmin, rmax, Nbin; logscale)
    r = abs(r)
    if rmin > r || r > rmax
        return nothing
    end
    if logscale
        dr = (log10(rmax) - log10(rmin))/Nbin
        return ceil(Int, (log10(r) - log10(rmin))/dr)
    else
        dr = (rmax - rmin)/Nbin
        return ceil(Int, (r - rmin)/dr)
    end
end

"""
    fold_histogram(hist)

    # Input
    - hist::Array{Float, 5}: A histogram of triplet counts for r12, r23, r31

    # Output
    - hist_red::Array{Float, 6}: A reduced histogram where different permutations of r12, r23, r31 
      are combined. Has 4 columns. First 5 for bin indexes (r12, r23, r31, π1, π2), the 6th for the total count.

    The condition imposed on the bins is r12 < r23 < r31.
"""
function fold_histogram(hist)
    Nrp = size(hist)[1]
    Nπ = size(hist)[4]
    Ntot = floor(Int, Nrp/6*(Nrp^2 + 3*Nrp + 2))*Nπ*Nπ
    hist_red = zeros(6, Ntot)
    count = 1
    for i in 1:Nrp, j in i:Nrp, k in j:Nrp, l in 1:Nπ, m in 1:Nπ
        hist_red[1,count] = i
        hist_red[2,count] = j
        hist_red[3,count] = k
        hist_red[4,count] = l
        hist_red[5,count] = m
        if i == j == k
            hist_red[6,count] = hist[i,j,k,l,m]
        elseif i == j
            hist_red[6,count] = hist[i,j,k,l,m] + hist[i,k,j,l,m] + hist[k,i,j,l,m]
        elseif j == k
            hist_red[6,count] = hist[i,j,k,l,m] + hist[k,i,j,l,m] + hist[k,j,i,l,m]
        elseif k == i
            hist_red[6,count] = hist[i,j,k,l,m] + hist[i,k,j,l,m] + hist[j,i,k,l,m]
        else
            hist_red[6,count] = hist[i,j,k,l,m] + hist[i,k,j,l,m] + hist[j,i,k,l,m] + hist[j,k,i,l,m] + hist[k,i,j,l,m] + hist[k,j,i,l,m]
        end
        count += 1
    end
    return hist_red
end