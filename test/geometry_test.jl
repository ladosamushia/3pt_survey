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