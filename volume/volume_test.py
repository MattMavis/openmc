#!/usr/bin/env python3

import openmc
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np
import datetime
from string import Formatter

mats = openmc.Materials()



#13%H2O+87%WC
H2OWC = openmc.Material(material_id=1,name='H2OWC')
H2OWC.set_density('g/cm3', 13.7281)
H2OWC.add_nuclide('H1',   0.0899, percent_type='ao')
H2OWC.add_nuclide('H2',   0.00001, percent_type='ao')
H2OWC.add_nuclide('O16',  0.0450, percent_type='ao')
H2OWC.add_nuclide('W182', 0.1142, percent_type='ao')
H2OWC.add_nuclide('W183', 0.0623, percent_type='ao')
H2OWC.add_nuclide('W184', 0.1324, percent_type='ao')
H2OWC.add_nuclide('W186', 0.1229, percent_type='ao')
H2OWC.add_nuclide('C12',  0.4327, percent_type='ao')
mats.append(H2OWC)

#CuCrZr@5.00e-02/Eurofer@9.00e-01/Helium@5.00e-02/Void@0.00e+00/Void@0.00e+00/Void@0.00e+00:
CuCrZr = openmc.Material(material_id=2,name='CuCrZr')
CuCrZr.set_density('g/cm3', 7.528)
CuCrZr.add_nuclide('Fe54',  4.82666e-02, percent_type='ao')
CuCrZr.add_nuclide('Fe56',  7.63944e-01, percent_type='ao')
CuCrZr.add_nuclide('Fe57',  1.74758e-02, percent_type='ao')
CuCrZr.add_nuclide('Fe58',  2.49655e-03, percent_type='ao')
CuCrZr.add_nuclide('Cr50',  4.16488e-03, percent_type='ao')
CuCrZr.add_nuclide('Cr52',  8.11668e-02, percent_type='ao')
CuCrZr.add_nuclide('Cr53',  9.20149e-03, percent_type='ao')
CuCrZr.add_nuclide('Cr54',  2.32459e-03, percent_type='ao')
CuCrZr.add_nuclide('C12',   5.26520e-03, percent_type='ao')
CuCrZr.add_nuclide('Mn55',  5.75567e-03, percent_type='ao')
CuCrZr.add_nuclide('V50',   6.46586e-06, percent_type='ao')
CuCrZr.add_nuclide('V51',   2.57988e-03, percent_type='ao')
CuCrZr.add_nuclide('N14',   1.68695e-03, percent_type='ao')
CuCrZr.add_nuclide('N15',   6.19691e-06, percent_type='ao')
CuCrZr.add_nuclide('O16',   3.29352e-04, percent_type='ao')
CuCrZr.add_nuclide('O17',   1.25455e-07, percent_type='ao')
CuCrZr.add_nuclide('Ta181', 2.62124e-04, percent_type='ao')
CuCrZr.add_nuclide('W182',  9.04794e-04, percent_type='ao')
CuCrZr.add_nuclide('W183',  4.92522e-04, percent_type='ao')
CuCrZr.add_nuclide('W184',  1.05737e-03, percent_type='ao')
CuCrZr.add_nuclide('W186',  9.85044e-04, percent_type='ao')
CuCrZr.add_nuclide('Cu63',  3.56780e-02, percent_type='ao')
CuCrZr.add_nuclide('Cu65',  1.58798e-02, percent_type='ao')
CuCrZr.add_nuclide('Zr90',  2.80384e-05, percent_type='ao')
CuCrZr.add_nuclide('Zr91',  6.09768e-06, percent_type='ao')
CuCrZr.add_nuclide('Zr92',  9.30985e-06, percent_type='ao')
CuCrZr.add_nuclide('Zr94',  9.47318e-06, percent_type='ao')
CuCrZr.add_nuclide('Zr96',  1.52442e-06, percent_type='ao')
CuCrZr.add_nuclide('He3',   2.10019e-11, percent_type='ao')
CuCrZr.add_nuclide('He4',   1.50014e-05, percent_type='ao')
mats.append(CuCrZr)

#Beryllium@5.50e-01/Eurofer@1.00e-01/Helium@2.00e-01/Lithium-Silicate@1.50e-01
Beryllium = openmc.Material(material_id=3,name='Beryllium')
Beryllium.set_density('g/cm3', 2.157)
Beryllium.add_nuclide('He4',   5.4569e-05, percent_type='ao')
Beryllium.add_nuclide('Li6',   6.2150e-03, percent_type='ao')
Beryllium.add_nuclide('Li7',   7.6652e-02, percent_type='ao')
Beryllium.add_nuclide('Be9',   7.2948e-01, percent_type='ao')
Beryllium.add_nuclide('C12',   4.1948e-04, percent_type='ao')
Beryllium.add_nuclide('N14',   1.0752e-04, percent_type='ao')
Beryllium.add_nuclide('N15',   3.9497e-07, percent_type='ao')
Beryllium.add_nuclide('O16',   7.7649e-02, percent_type='ao')
Beryllium.add_nuclide('Si28',  1.7924e-02, percent_type='ao')
Beryllium.add_nuclide('Si29',  9.0759e-04, percent_type='ao')
Beryllium.add_nuclide('Si30',  6.0247e-04, percent_type='ao')
Beryllium.add_nuclide('V50',   1.9781e-04, percent_type='ao')
Beryllium.add_nuclide('Cr50',  3.7893e-04, percent_type='ao')
Beryllium.add_nuclide('Cr52',  7.3072e-03, percent_type='ao')
Beryllium.add_nuclide('Cr53',  8.2858e-04, percent_type='ao')
Beryllium.add_nuclide('Cr54',  2.0625e-04, percent_type='ao')
Beryllium.add_nuclide('Mn55',  3.6689E-04, percent_type='ao')
Beryllium.add_nuclide('Fe54',  4.7397e-03, percent_type='ao')
Beryllium.add_nuclide('Fe56',  7.3683e-02, percent_type='ao')
Beryllium.add_nuclide('Fe57',  1.6870e-03, percent_type='ao')
Beryllium.add_nuclide('Fe58',  2.2493e-04, percent_type='ao')
Beryllium.add_nuclide('Ta181', 3.3414e-05, percent_type='ao')
Beryllium.add_nuclide('W182',  7.9287e-05, percent_type='ao')
Beryllium.add_nuclide('W183',  4.3111e-05, percent_type='ao')
Beryllium.add_nuclide('W184',  9.2552e-05, percent_type='ao')
Beryllium.add_nuclide('W186',  8.6221e-05, percent_type='ao')
mats.append(Beryllium)

#70%EUROFER 30%Water
eurofer70_H2030 = openmc.Material(material_id=4,name='eurofer70_H2030')
eurofer70_H2030.set_density('g/cm3', 5.758)
eurofer70_H2030.add_nuclide('H1',    2.2373e-01, percent_type='ao')
eurofer70_H2030.add_nuclide('H2',    3.3565e-05, percent_type='ao')
eurofer70_H2030.add_nuclide('C12',   3.0751e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('N14',   7.8819E-03, percent_type='ao')
eurofer70_H2030.add_nuclide('N15',   2.8954e-06, percent_type='ao')
eurofer70_H2030.add_nuclide('O16',   1.1290e-01, percent_type='ao')
eurofer70_H2030.add_nuclide('V50',   1.4501e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Cr50',  2.7778e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Cr52',  5.3567e-02, percent_type='ao')
eurofer70_H2030.add_nuclide('Cr53',  6.0741e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Cr54',  1.5120e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Mn55',  2.6896e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Fe54',  3.4745e-02, percent_type='ao')
eurofer70_H2030.add_nuclide('Fe56',  5.4014e-01, percent_type='ao')
eurofer70_H2030.add_nuclide('Fe57',  1.2367e-02, percent_type='ao')
eurofer70_H2030.add_nuclide('Fe58',  1.6489e-03, percent_type='ao')
eurofer70_H2030.add_nuclide('Ta181', 2.4494e-04, percent_type='ao')
eurofer70_H2030.add_nuclide('W182',  5.8123e-04, percent_type='ao')
eurofer70_H2030.add_nuclide('W183',  3.1603e-04, percent_type='ao')
eurofer70_H2030.add_nuclide('W184',  6.7847e-04, percent_type='ao')
eurofer70_H2030.add_nuclide('W186',  6.3206e-04, percent_type='ao')
mats.append(eurofer70_H2030)

#NB3SN BASED STRAND
# MASS DENSITY [G/CC] - 8.94
# VOLUME FRACTION [%] - 13.86
# T.A.D. = 1.08843E-002
# EFF.DENSITY = 1.23908E+000
NB3SN = openmc.Material(material_id=5,name='NB3SN')
NB3SN.set_density('g/cm3', 8.94)
NB3SN.add_nuclide('Cu63',  6.48166e-03, percent_type='ao')
NB3SN.add_nuclide('Cu65',  2.88896e-03, percent_type='ao')
NB3SN.add_nuclide('Nb93',  9.63801e-04, percent_type='ao')
NB3SN.add_nuclide('Sn112', 2.43891e-06, percent_type='ao')
NB3SN.add_nuclide('Sn114', 1.65946e-06, percent_type='ao')
NB3SN.add_nuclide('Sn115', 8.54876e-07, percent_type='ao')
NB3SN.add_nuclide('Sn116', 3.65585e-05, percent_type='ao')
NB3SN.add_nuclide('Sn117', 1.93101e-05, percent_type='ao')
NB3SN.add_nuclide('Sn118', 6.08973e-05, percent_type='ao')
NB3SN.add_nuclide('Sn119', 2.15982e-05, percent_type='ao')
NB3SN.add_nuclide('Sn120', 8.19172e-05, percent_type='ao')
NB3SN.add_nuclide('Sn122', 1.16414e-05, percent_type='ao')
NB3SN.add_nuclide('Sn124', 1.45580e-05, percent_type='ao')
NB3SN.add_nuclide('Ta181', 1.23714e-04, percent_type='ao')
NB3SN.add_nuclide('Ti46',  2.57138e-06, percent_type='ao')
NB3SN.add_nuclide('Ti47',  2.31892e-06, percent_type='ao')
NB3SN.add_nuclide('Ti48',  2.29773e-05, percent_type='ao')
NB3SN.add_nuclide('Ti49',  1.68620e-06, percent_type='ao')
NB3SN.add_nuclide('Ti50',  1.61452e-06, percent_type='ao')
NB3SN.add_nuclide('Cr50',  6.23551e-06, percent_type='ao')
NB3SN.add_nuclide('Cr52',  1.20246e-04, percent_type='ao')
NB3SN.add_nuclide('Cr53',  1.36349e-05, percent_type='ao')
NB3SN.add_nuclide('Cr54',  3.39401e-06, percent_type='ao')
mats.append(NB3SN)

mats.export_to_xml()

#define all the surfaces
central_sol_surface = openmc.ZCylinder(r=50, boundary_type='transmission') #50 Centre column radius (cm)

central_shield_outer_surface = openmc.ZCylinder(r=90, boundary_type='transmission') #40 Centre column n-shield thickness (cm)

first_wall_centre = openmc.ZCylinder(r=92, boundary_type='vacuum') #First wall centre collumn 2cm thick

first_wall_inner_surface = openmc.Sphere(r=420.2, boundary_type='vacuum') #First wall inner sphere

first_wall_outer_surface = openmc.Sphere(r=422.2, boundary_type='transmission') #First wall outer sphere

breeder_blanket_outer_surface = openmc.Sphere(r=492.2, boundary_type='transmission') #Breeder blanket 70cm thick

outer_sphere = openmc.Sphere(r=522.2, boundary_type='vacuum') #Outer sphere

surface_100 = openmc.ZPlane(z0=-5.0, boundary_type='transmission')
surface_101 = openmc.ZPlane(z0=5.0, boundary_type='transmission')
surface_200 = openmc.XPlane(x0=0.0, boundary_type='vacuum')
surface_201 = openmc.YPlane(y0=0.0, boundary_type='vacuum')

#define the cells
#Central Solenoid
central_sol_region_1 = -central_sol_surface & +surface_101 & -outer_sphere & +surface_200 & -surface_201
central_sol_cell_1 = openmc.Cell(region=central_sol_region_1)
central_sol_cell_1.fill = NB3SN

central_sol_region_2 = -central_sol_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
central_sol_cell_2 = openmc.Cell(region=central_sol_region_2)
central_sol_cell_2.fill = NB3SN

central_sol_region_3 = -central_sol_surface & -surface_100 & -outer_sphere & +surface_200 & -surface_201
central_sol_cell_3 = openmc.Cell(region=central_sol_region_3)
central_sol_cell_3.fill = NB3SN

#Central Solenoid Shield
central_sol_shield_region_1 = -central_shield_outer_surface & +central_sol_surface & +surface_101 & -outer_sphere & +surface_200 & -surface_201
central_sol_shield_cell_1 = openmc.Cell(region=central_sol_shield_region_1)
central_sol_shield_cell_1.fill = H2OWC

central_sol_shield_region_2 = -central_shield_outer_surface & +central_sol_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
central_sol_shield_cell_2 = openmc.Cell(region=central_sol_shield_region_2)
central_sol_shield_cell_2.fill = H2OWC

central_sol_shield_region_3 = -central_shield_outer_surface & +central_sol_surface & -surface_100 & -outer_sphere & +surface_200 & -surface_201
central_sol_shield_cell_3 = openmc.Cell(region=central_sol_shield_region_3)
central_sol_shield_cell_3.fill = H2OWC

#Central Solenoid First Wall
first_wall_centre_region_1 = -first_wall_centre & +central_shield_outer_surface & +surface_101 & -outer_sphere & +surface_200 & -surface_201
first_wall_centre_cell_1 = openmc.Cell(region=first_wall_centre_region_1)
first_wall_centre_cell_1.fill = CuCrZr

first_wall_centre_region_2 = -first_wall_centre & +central_shield_outer_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
first_wall_centre_cell_2 = openmc.Cell(region=first_wall_centre_region_2)
first_wall_centre_cell_2.fill = CuCrZr

first_wall_centre_region_3 = -first_wall_centre & +central_shield_outer_surface & -surface_100 & -outer_sphere & +surface_200 & -surface_201
first_wall_centre_cell_3 = openmc.Cell(region=first_wall_centre_region_3)
first_wall_centre_cell_3.fill = CuCrZr

#Outer First Wall
first_wall_outer_region_1 = -first_wall_outer_surface & +first_wall_inner_surface & +first_wall_centre & +surface_101 & +surface_200 & -surface_201
first_wall_outer_cell_1 = openmc.Cell(region=first_wall_outer_region_1)
first_wall_outer_cell_1.fill = CuCrZr

first_wall_outer_region_2 = -first_wall_outer_surface & +first_wall_inner_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
first_wall_outer_cell_2 = openmc.Cell(region=first_wall_outer_region_2)
first_wall_outer_cell_2.fill = CuCrZr 

first_wall_outer_region_3 = -first_wall_outer_surface & +first_wall_inner_surface & +first_wall_centre & -surface_100 & +surface_200 & -surface_201
first_wall_outer_cell_3 = openmc.Cell(region=first_wall_outer_region_3)
first_wall_outer_cell_3.fill = CuCrZr 

#Breeder Blanket
breeder_blanket_region_1 = -breeder_blanket_outer_surface & +first_wall_outer_surface & +first_wall_centre & +surface_101 & +surface_200 & -surface_201
breeder_blanket_cell_1 = openmc.Cell(region=breeder_blanket_region_1)
breeder_blanket_cell_1.fill = Beryllium

breeder_blanket_region_2 = -breeder_blanket_outer_surface & +first_wall_outer_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
breeder_blanket_cell_2 = openmc.Cell(region=breeder_blanket_region_2)
breeder_blanket_cell_2.fill = Beryllium

breeder_blanket_region_3 = -breeder_blanket_outer_surface & +first_wall_outer_surface & +first_wall_centre & -surface_100 & +surface_200 & -surface_201
breeder_blanket_cell_3 = openmc.Cell(region=breeder_blanket_region_3)
breeder_blanket_cell_3.fill = Beryllium

#Cyrogen
cyrogen_region_1 = -outer_sphere & +breeder_blanket_outer_surface & +first_wall_centre & +surface_101 & +surface_200 & -surface_201
cyrogen_cell_1 = openmc.Cell(region=cyrogen_region_1)
cyrogen_cell_1.fill = eurofer70_H2030

cyrogen_region_2 = -outer_sphere & +breeder_blanket_outer_surface & -surface_101 & +surface_100 & +surface_200 & -surface_201
cyrogen_cell_2 = openmc.Cell(region=cyrogen_region_2)
cyrogen_cell_2.fill = eurofer70_H2030

cyrogen_region_3 = -outer_sphere & +breeder_blanket_outer_surface & +first_wall_centre & -surface_100 & +surface_200 & -surface_201
cyrogen_cell_3 = openmc.Cell(region=cyrogen_region_3)
cyrogen_cell_3.fill = eurofer70_H2030

#Voids
inner_void_region = -first_wall_inner_surface & +first_wall_centre & +surface_200 & -surface_201
inner_void_cell = openmc.Cell(region=inner_void_region)

outer_void_region = +outer_sphere | -surface_200 | +surface_201
outer_void_cell = openmc.Cell(region=outer_void_region)

universe = openmc.Universe(0,name='Main_Universe', cells=[central_sol_cell_1, central_sol_cell_2, central_sol_cell_3, central_sol_shield_cell_1, central_sol_shield_cell_2, central_sol_shield_cell_3, first_wall_centre_cell_1, first_wall_centre_cell_2, first_wall_centre_cell_3, first_wall_outer_cell_1, first_wall_outer_cell_2, first_wall_outer_cell_3, breeder_blanket_cell_1, breeder_blanket_cell_2, breeder_blanket_cell_3, cyrogen_cell_1, cyrogen_cell_2, cyrogen_cell_3, inner_void_cell, outer_void_cell])
geom = openmc.Geometry(universe)
geom.export_to_xml()
print ("geo export complete")
print ("defining volumes to calculate")
#Settings
###################################
# Write comments for this section #
###################################
settings = openmc.Settings()
lower_left_x_int = 0
lower_left_y_int = -524
lower_left_z_int = -524
upper_right_x_limit = 524 
upper_right_y_limit = 0
upper_right_z_limit = 524
pitch = 262
vol_calc_array = []
resultCoord = []
j = 0
lower_left_x = lower_left_x_int
while lower_left_x<=upper_right_x_limit:
    lower_left_y = lower_left_y_int
    while lower_left_y<=upper_right_y_limit:
        lower_left_z = lower_left_z_int
        while lower_left_z<=upper_right_z_limit:
            lower_left = (lower_left_x, lower_left_y, lower_left_z)
            if lower_left_x >= upper_right_x_limit:
                break
            if lower_left_y >= upper_right_y_limit:
                break
            if lower_left_z >= upper_right_z_limit:
                break
            
            upper_right_x = lower_left_x + pitch
            upper_right_y = lower_left_y + pitch
            upper_right_z = lower_left_z + pitch
            upper_right = (upper_right_x, upper_right_y, upper_right_z)
            mid_coord = (lower_left_x + (pitch/2),lower_left_y + (pitch/2),lower_left_z + (pitch/2))
            #print(lower_left)
            #print(upper_right)
            resultCoord.append(mid_coord)
            vol_calc = openmc.VolumeCalculation([H2OWC, CuCrZr, Beryllium, eurofer70_H2030, NB3SN], 1000000, lower_left, upper_right)
            #H2OWC, CuCrZr, Beryllium, eurofer70_H2030, NB3SN
            #vol_calc = openmc.VolumeCalculation([central_sol_cell_1, central_sol_cell_2, central_sol_cell_3, central_sol_shield_cell_1, central_sol_shield_cell_2, central_sol_shield_cell_3, first_wall_centre_cell_1, first_wall_centre_cell_2, first_wall_centre_cell_3, first_wall_outer_cell_1, first_wall_outer_cell_2, first_wall_outer_cell_3, breeder_blanket_cell_1, breeder_blanket_cell_2, breeder_blanket_cell_3, cyrogen_cell_1, cyrogen_cell_2, cyrogen_cell_3, inner_void_cell, outer_void_cell], 1000000, lower_left, upper_right)
            vol_calc_array.append(vol_calc)
            j += 1
            lower_left_z += pitch
        lower_left_y += pitch
    lower_left_x += pitch
#output bounds#
xbound_value = lower_left_x_int
ybound_value = lower_left_y_int
zbound_value = lower_left_z_int
bounds_x = [xbound_value]
bounds_y = [ybound_value]
bounds_z = [zbound_value]
while xbound_value < upper_right_x_limit:
    xbound_value += pitch
    bounds_x.append(xbound_value)

while ybound_value < upper_right_y_limit:
    ybound_value += pitch
    bounds_y.append(ybound_value)

while zbound_value < upper_right_z_limit:
    zbound_value += pitch
    bounds_z.append(zbound_value)

print(j)
settings.volume_calculations = vol_calc_array
print("defining source parameters")
#Source
batches = 2
settings.batches = batches
settings.inactive = 0
particles = 100000
settings.particles = particles
settings.run_mode = 'fixed source'
histories = (settings.particles*batches) - (settings.inactive*batches)
source = openmc.Source()
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])
source.space = openmc.stats.Point((266,-266,0))
settings.source = source
xml_source = source.to_xml_element()
settings.export_to_xml()
print("settings exported")
print("defining tally mesh")
#Mesh
mesh_upper_right = (upper_right_x_limit,upper_right_y_limit,upper_right_z_limit)
mesh_lower_left = (lower_left_x_int,lower_left_y_int,lower_left_z_int)
mesh = openmc.RegularMesh()
xdimension = len(bounds_x) -1
#print(int(xdimension))
ydimension = len(bounds_y) -1
#print(ydimension)
zdimension = len(bounds_z) -1
#print(zdimension)
mesh.dimension = [xdimension, ydimension, zdimension]
mesh.width = [pitch,pitch,pitch]
mesh.lower_left = mesh_lower_left
#mesh.upper_right = mesh_upper_right
#Tallies
tallies = openmc.Tallies()
#mesh filter
mesh_filter = openmc.MeshFilter(mesh)
mesh_height = pitch
mesh_width = pitch
mesh_tally = openmc.Tally(1,name='tallies_on_mesh')
energy_bins = openmc.mgxs.GROUP_STRUCTURES['VITAMIN-J-175']
energy_filter = openmc.EnergyFilter(energy_bins)
mesh_tally.filters = [mesh_filter,energy_filter]
tallies_to_plot = 'flux'
mesh_tally.scores = [tallies_to_plot]
tallies.append(mesh_tally)
tallies.export_to_xml()
print("Tallies exported")
print ("defining plot")
# makes the 3d "cube" style geometry 
vox_plot = openmc.Plot()
#plt.show(universe.plot(width=(2000,2000),basis='xz'))
#plt.show(universe.plot(width=(2000,2000),basis='xy'))
#plt.show(universe.plot(width=(2000,2000),basis='yz'))
vox_plot.type = 'voxel'
vox_plot.width = (1500., 1500., 1500.)
vox_plot.pixels = (300, 300, 300)
vox_plot.filename = 'plot_3d_vol_test'
vox_plot.color_by = 'material'
plots = openmc.Plots([vox_plot])
plots.export_to_xml()
print("plot exported")
#Run OpenMC
#print("Running OpenMC geometry plot")
#openmc.plot_geometry(openmc_exec='/home/mmavis/openmc/bin/openmc')
print("running OpenMC")
openmc.run(openmc_exec='/home/mmavis/openmc/bin/openmc')
#print("calculating volumes")
openmc.calculate_volumes(openmc_exec='/home/mmavis/openmc/bin/openmc')

#Following is the code for reading in the volume_*.h5 files and outputtin the results to the Output.txt file
numResults = j      #set the amount of results generated. j is defined in the code for creating the volume calculation settings
results = []        #initialise results array
#print(numResults)
i = 1

#Read in the volume resutls from files and add them to the results array
while i <= numResults:
    vol_calc.load_results('volume_' + str(i) +'.h5')
    #print(vol_calc.volumes)
    results.append(vol_calc.volumes)
    i+=1
#print(results)


####################
# File Output Code #
####################
#Turn into subroutine#
k = 0
l = 0
nint = [len(bounds_x),len(bounds_y),len(bounds_z)]
outputVolPer = []
#Calculate the total volume of the sampled areas.
total_volume = (pitch**3)                                           
coord = 'X Cord','Y Cord','Z Cord'                     
#Load the first result.
result = results[k]                                                 
#Find out the amount of values in the result.
numVolumes = len(result)                                            
#Break the results into the 2 seperate arrays.
mat, vol = zip(*result.items())                                     
#Open the Output.txt file and clear it if it already exists. If not create it.
open('posmat.txt', 'w').close()                                     
#Open Output.txt in append mode.
with open('posmat.txt', 'a+') as file_handler:
    #Output the amount of bounds in each direction
    for item in nint:
        file_handler.write('{0:>12d}'.format(item))
    file_handler.write("\n")
    #Output X bounds
    file_handler.write("xbounds ")
    for item in bounds_x:
        file_handler.write('{0:10.2f}'.format(item))
    file_handler.write("\n")
    #Output Y bounds
    file_handler.write("ybounds ")
    for item in bounds_y:
        file_handler.write('{0:10.2f}'.format(item))
    file_handler.write("\n")
    #Output Z bounds
    file_handler.write("zbounds ")
    for item in bounds_z:
        file_handler.write('{0:10.2f}'.format(item))
    file_handler.write("\n")
    #Iterate through the values in coord and write them to Output.txt, each value is given 12 characters and is justified to the right
    f=0
    for item in coord:
        fmt=['{0:>12s}','{0:>13s}','{0:>13s}']                                            
        file_handler.write(fmt[f].format(item))
        f += 1
    #Iterate through the values in mat and write them to Output.txt, each value is given 8 integers and is justified right.
    for item in mat:                                                
        file_handler.write(str(' ')+'{0:8d}'.format(item))
    #Write void to the end of the line with 8 characters and is justified right, then moves to a new line in preperation for the co-ordinates and volumes.
    file_handler.write(str(' ')+'{0:>8s}'.format(str('Void'))+"\n")
#Iterate through all the results
for k in range(numResults):
    #Reset l
    l = 0
    #Load result in array position k
    result = results[k]
    #Break the results into the 2 seperate arrays.
    mat, vol = zip(*result.items())
    #Find out the amount of values in the result.
    numVolumes = len(result)
    #Calculate the percentage of each material in the sampled area. (results range from values 0-1).
    outputVolPer =[(vol[l].n/(total_volume)) for l in range(numVolumes)]
    #Print statement for debug
    [print('Material ' +str(mat[l]) + ': ' + str('{0:.4f}'.format(outputVolPer[l]*100)) +'%') for l in range(numVolumes)]
    #Calculate the amount of void present in the sampled area.
    void = 1 - sum(outputVolPer)
    #Print statement for debug.
    print('Void: ' + str('{0:.4f}'.format(void*100))+ '%')
    #Add void value to the end of the volume array.
    outputVolPer.append(void)
    #Set the coordinate the volumes were taken at. Defined on line 294.
    resultCoordOut = resultCoord[k]
    #Open Output.txt in append mode.
    with open('posmat.txt', 'a+') as file_handler:
        #Iterate through the values in resultCoordOut and write them to Output.txt, each value is in scientific form with 4 decimal places and is justified to the right.
        f=0
        for item in resultCoordOut:
            fmt=["{0:12.4e}","{0:13.4e}","{0:13.4e}"]
            file_handler.write(fmt[f].format(item))
            f += 1
        #Iterate through the values in outputVolPer and write them to Output.txt, each value is a float with 8 places with 6 of them decimal places and is justified to the right.
        for item in outputVolPer:
            file_handler.write(str(" ")+"{0:8.6f}".format(item))
        #Start new line for next result.
        file_handler.write('\n')
    #Close the file
    file_handler.close()
#Delete the volume_*.h5 files to clean up the directory.
os.system('rm -r volume_*.h5')

#Not working currently
sp = openmc.StatePoint('statepoint.'+str(batches)+'.h5')
flux_tally = sp.get_tally(scores=[tallies_to_plot])
tally_results = flux_tally.get_pandas_dataframe()
tally_number = 1 #figure out a way to change this to array of ids
print(tally_results)
n=1
for n in range(len(bounds_x)):
    tally_results[('mesh 1','x')].replace(n,(lower_left_x_int+(pitch*(n-1)+(pitch/2))),inplace=True)

n=1
for n in range(len(bounds_y)):
    tally_results[('mesh 1','y')].replace(n,(lower_left_y_int+(pitch*(n-1)+(pitch/2))),inplace=True)                

n=1
for n in range(len(bounds_z)):
    tally_results[('mesh 1','z')].replace(n,(lower_left_z_int+(pitch*(n-1)+(pitch/2))),inplace=True)

tally_results.loc[:,'X'] = tally_results[('mesh 1','x')]
tally_results.loc[:,'Y'] = tally_results[('mesh 1','y')]
tally_results.loc[:,'Z'] = tally_results[('mesh 1','z')]
tally_results = tally_results.drop("energy high [eV]",axis=1)
tally_results = tally_results.drop("nuclide",axis=1)
tally_results = tally_results.drop("score",axis=1)
tally_results = tally_results.drop("mesh 1",axis=1)
tally_results = tally_results[['energy low [eV]','X','Y','Z','mean','std. dev.',]]
tally_results.head()
tally_results.rename(columns={'energy low [eV]': 'Energy', 'mean':'Results', 'std. dev.': 'Rel Error'}, inplace=True)
tally_results.sort_values(by=['Energy','X','Y','Z'],ascending=[True,True,True,True],inplace=True)
tally_results.reset_index(inplace=True,drop=True)

#print(tally_results)
#with open('tally_results.txt', 'w') as tally_file:
#    tally_file.write(tally_results.to_string(index=False))
#    tally_file.close()
#tally_results = tally_results.to_dict()
#print(tally_results)
#print(energy_bins)
energy_bins_np = np.array([energy_bins])
tally_results = tally_results.to_numpy()
#print(energy_bins_np)
tally_titles = np.array([['Energy','X','Y','Z','Results','Rel Error']])
print(histories)

open('meshtal.msht', 'w').close() 
with open('meshtal.msht','a+') as tally_file:
    today = datetime.datetime.now()
    tally_file.write("OpenMC   version 0.11.0  "+today.strftime('%d/%m/%Y')+ " probid =  "+today.strftime('%d/%m/%Y %H:%M:%S'+"\n"))
    tally_file.write(" Cube"+"\n")
    tally_file.write(" Number of histories used for normalizing tallies ="+"{0:17.2f}".format(histories)+"\n"+"\n")
    tally_file.write(" Mesh Tally Number"+"{0:10d}".format(tally_number)+"\n"+" neutron  mesh tally." +"\n"+"\n")
    tally_file.write(" Tally bin boundaries:"+"\n")
    #Output X bounds
    tally_file.write("    X direction:")
    for item in bounds_x:
        tally_file.write('{0:10.2f}'.format(item))
    tally_file.write("\n")
    #Output Y bounds
    tally_file.write("    Y direction:")
    for item in bounds_y:
        tally_file.write('{0:10.2f}'.format(item))
    tally_file.write("\n")
    #Output Z bounds
    tally_file.write("    Z direction:")
    for item in bounds_z:
        tally_file.write('{0:10.2f}'.format(item))
    tally_file.write("\n")
    tally_file.write("    Energy bin boundaries: ")
    np.savetxt(tally_file,(energy_bins_np),fmt="%8.2e",newline='\n')
    tally_file.write("\n")
    np.savetxt(tally_file,tally_titles,fmt="%11s %9s %9s %9s %11s %11s ",delimiter=',',newline='\n')
    np.savetxt(tally_file,tally_results,fmt="%11.3e %9.3f %9.3f %9.3f %11.5e %11.5e",delimiter=' ',newline='\n')
print(xml_source)
# this converts the h5 file to a vti
#os.system('/home/mmavis/openmc/scripts/openmc-voxel-to-silovtk plot_3d_vol_test.h5 -o plot_3d_vol_test.vti')
