#
#
#

"""
    par_perp_distance(xyz1, xyz2)

    # Input
    - xyz1::Array{Float, 1}: x, y, z coordinates of point 1
    - xyz2::Array{Float, 1}: x, y, z coordinates of point 2

    # Output
    - ::Float: parallel distance with respect to the line of sight
    - ::Float: Perpendicular distance with respect to the line of sight
"""
function par_perp_distance(xyz1, xyz2)
    # Vector connecting 1 to 2
    dxyz = xyz1 .- xyz2
    # Vector pointing in between 1 and 2
    mxyz = (xyz1 .+ xyz2)./2
    # Project onto mxyz
    xyz_par = mxyz.*dxyz/sqrt(sum(mxyz.^2))
    # Perpendicular to mxyz
    xyz_per = dxyz .- xyz_par
    return sqrt(sum(xyz_par.^2)), sqrt(sum(xyz_per.^2))
end

"""
    radecz_to_xyz(ra, dec, z)

    # Input
    - ra::
"""
function radecz_to_xyz(ra, dec, z)
    Ngal = length(z)
    xyz = zeros(3, Ngal)
    for i in 1:Ngal
        x[1,i] = 
        y[2,i] = 
        z[3,i] = 
    return xyz
end