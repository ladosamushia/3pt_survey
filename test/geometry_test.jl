using Test

include("../src/geometry.jl")

xyz1 = [5, 0, 100]
xyz2 = [-5, 0, 100]
dpar, dper = par_perp_distance(xyz1, xyz2)
@test dpar == 0
@test dper == 10
dpar, dper = par_perp_distance(xyz2, xyz1)
@test dpar == 0
@test dper == 10
xyz1 = [500, 1, 1]
xyz2 = [500, -1, -1]
dpar, dper = par_perp_distance(xyz1, xyz2)
@test dpar == 0
@test isapprox(dper, 2*sqrt(2))
xyz1 = [0, 0, 102]
xyz2 = [0, 0, 98]
dpar, dper = par_perp_distance(xyz1, xyz2)
@test dpar == 4
@test dper == 0
xyz1 = [5*cos(1), 5*sin(1), 100]
xyz2 = [-5*cos(1), -5*sin(1), 100]
dpar, dper = par_perp_distance(xyz1, xyz2)
@test dpar == 0
@test dper == 10

ra = [0, 2*π*0.136, 2*π*0.84]
dec = [0, π/2, π/3]
z = [1, 0.5, 1.5]
xyz = radecz_to_xyz(ra, dec, z)
@test xyz[1,:] .≊ [3433.4, 0, 0]
@test xyz[2,:] .≊ []
@test xyz[3,:] .≊ []