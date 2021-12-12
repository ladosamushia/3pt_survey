#
#
#

using QuadGK

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
    - ra::Array{Float, 1}: right ascention in radians
    - dec::Array{Float, 1}: declanation in radians
    - z::Array{Float, 1}: redshift

    # Output
    - xyz::Array{Float, 3}: Cartesian x, y, z, coordinates

    The AbacusSummit cosmology is hardcoded into this.
"""
function radecz_to_xyz(ra, dec, z)
    # Ωm in AbacusSummit
    Ωm = 0.31377
    h = 0.6736
    Ngal = length(z)
    r = zeros(Ngal)
    xyz = zeros(3, Ngal)
    for i in 1:Ngal
        dist, _ = quadgk(x -> 2997.92458/h/sqrt(Ωm*(1 + x)^3 + 1 - Ωm), 0, z[i], rtol=1e-5)
        r[i] = dist
    end
    xyz[1,:] = r.*cos.(ra).*cos.(dec)
    xyz[2,:] = r.*sin.(ra).*cos.(dec)
    xyz[3,:] = r.*sin.(dec)
    return xyz
end