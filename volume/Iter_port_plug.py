#!/usr/bin/env python3

import openmc
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import figure, cm
import os
import h5py
import numpy as np
import datetime
from string import Formatter
from sys import getsizeof

mats = openmc.Materials()
#Steel and water
Steel_Water = openmc.Material(material_id=1,name='Steel_Water')
Steel_Water.set_density('g/cm3', 6.536)
Steel_Water.add_element('H', 1.46e-01, percent_type='ao')
Steel_Water.add_element('B', 4.02e-05, percent_type='ao')
Steel_Water.add_element('C', 8.14e-04, percent_type='ao')
Steel_Water.add_element('N', 2.17e-03, percent_type='ao')
Steel_Water.add_element('O', 7.29e-02, percent_type='ao')
Steel_Water.add_element('Al', 8.04e-04, percent_type='ao')
Steel_Water.add_element('Si', 7.73e-03, percent_type='ao')
Steel_Water.add_element('P', 3.50E-04, percent_type='ao')
Steel_Water.add_element('S', 1.02e-04, percent_type='ao')
Steel_Water.add_element('K', 5.55e-06, percent_type='ao')
Steel_Water.add_element('Ti', 1.36e-03, percent_type='ao')
Steel_Water.add_element('V', 3.41e-05, percent_type='ao')
Steel_Water.add_element('Cr', 1.46e-01, percent_type='ao')
Steel_Water.add_element('Mn', 1.42e-02, percent_type='ao')
Steel_Water.add_element('Fe', 5.03e-01, percent_type='ao')
Steel_Water.add_element('Co', 3.68e-04, percent_type='ao')
Steel_Water.add_element('Ni', 9.06e-02, percent_type='ao')
Steel_Water.add_element('Cu', 2.05e-03, percent_type='ao')
Steel_Water.add_element('Zr', 9.52e-06, percent_type='ao')
Steel_Water.add_element('Nb', 9.52e-06, percent_type='ao')
Steel_Water.add_element('Mo', 1.13e-02, percent_type='ao')
Steel_Water.add_element('Sn', 7.31e-06, percent_type='ao')
Steel_Water.add_element('Ta', 2.40e-05, percent_type='ao')
Steel_Water.add_element('W', 2.36e-06, percent_type='ao')
Steel_Water.add_element('Pb', 1.68e-06, percent_type='ao')
Steel_Water.add_element('Bi', 1.66e-06, percent_type='ao')
mats.append(Steel_Water)

#Steel
Steel = openmc.Material(material_id=2,name='Steel')
Steel.set_density('g/cm3', 7.93)
Steel.add_element('B', 5.14e-05, percent_type='ao')
Steel.add_element('C', 1.04e-03, percent_type='ao')
Steel.add_element('N', 2.78e-03, percent_type='ao')
Steel.add_element('O', 6.95e-05, percent_type='ao')
Steel.add_element('Al', 1.03e-03, percent_type='ao')
Steel.add_element('Si', 9.89e-03, percent_type='ao')
Steel.add_element('P', 4.48e-04, percent_type='ao')
Steel.add_element('S', 1.30e-04, percent_type='ao')
Steel.add_element('K', 7.10e-06, percent_type='ao')
Steel.add_element('Ti', 1.74e-03, percent_type='ao')
Steel.add_element('V', 4.36e-05, percent_type='ao')
Steel.add_element('Cr', 1.87e-01, percent_type='ao')
Steel.add_element('Mn', 1.82e-02, percent_type='ao')
Steel.add_element('Fe', 6.44e-01, percent_type='ao')
Steel.add_element('Co', 4.71e-04, percent_type='ao')
Steel.add_element('Ni', 1.16e-01, percent_type='ao')
Steel.add_element('Cu', 2.62e-03, percent_type='ao')
Steel.add_element('Zr', 1.22e-05, percent_type='ao')
Steel.add_element('Nb', 5.98e-05, percent_type='ao')
Steel.add_element('Mo', 1.45e-02, percent_type='ao')
Steel.add_element('Sn', 9.36e-06, percent_type='ao')
Steel.add_element('Ta', 3.07e-05, percent_type='ao')
Steel.add_element('W', 3.02e-06, percent_type='ao')
Steel.add_element('Pb', 2.14e-06, percent_type='ao')
Steel.add_element('Bi', 2.13e-06, percent_type='ao')
mats.append(Steel)

mats.export_to_xml()

#plug hole
plug_hole = openmc.YCylinder(r=7.5)

#Steel & water cylinder
SW_cylinder = openmc.YCylinder(r=48)

#inner cylinder
inner_cylinder = openmc.YCylinder(r=50)

#outer cylinder
outer_cylinder = openmc.YCylinder(r=100, boundary_type='vacuum')

#Plug Front
plug_front = openmc.YPlane(y0=110)

#Plug Back
plug_back = openmc.YPlane(y0=660)

#Steel and Water back
SW_back = openmc.YPlane(y0=320)

#inner cylinder back
inner_cylinder_back = openmc.YPlane(y0=645)

#Boundary Front
boundary_front = openmc.YPlane(y0=0, boundary_type='vacuum')

#Boundary back
boundary_back = openmc.YPlane(y0=700, boundary_type='vacuum')


#Tally Cell Front
tally_cell_front = openmc.YPlane(y0=690)

#Tally Cell 1
tally_cell_1_cylinder = openmc.YCylinder(r=15)

#Tally Cell 2
tally_cell_2_cylinder = openmc.YCylinder(r=30) 

#Tally Cell 3
tally_cell_3_cylinder = openmc.YCylinder(r=45)

#Tally Cell 4
tally_cell_4_cylinder = openmc.YCylinder(r=60)

#Source back
source_back = openmc.YPlane(y0=10)



#cells
#Plug hole cell
plug_hole_region = -plug_hole & +plug_front & -SW_back
plug_hole_cell = openmc.Cell(region=plug_hole_region)

#Steel & water cylinder
SW_Cylinder_region = +plug_hole & -SW_cylinder & +plug_front & -SW_back
SW_cylinder_cell = openmc.Cell(region=SW_Cylinder_region)
SW_cylinder_cell.fill = Steel_Water


#Steel back plate
steel_back_plate_region = -SW_cylinder & +SW_back & +inner_cylinder_back & -plug_back
steel_back_plate_cell = openmc.Cell(region=steel_back_plate_region)
steel_back_plate_cell.fill = Steel

#Inner cylinder
inner_cylinder_region = -SW_cylinder & +SW_back & -inner_cylinder_back
inner_cylinder_cell = openmc.Cell(region=inner_cylinder_region)

#Void gap between inner and outer cylinders
gap_2cm_region = +SW_cylinder & -inner_cylinder & +plug_front & -plug_back
gap_2cm_cell = openmc.Cell(region=gap_2cm_region)

#Steel cylinder
steel_cylinder_region = +inner_cylinder & -outer_cylinder & +plug_front & -plug_back
steel_cylinder_cell = openmc.Cell(region=steel_cylinder_region)
steel_cylinder_cell.fill = Steel

#Source cell
source_region = -outer_cylinder & +boundary_front & -source_back
source_cell = openmc.Cell(region=source_region)

#Front void gap
front_void_gap_region = -outer_cylinder & +source_back & -plug_front
front_void_gap_cell = openmc.Cell(region=front_void_gap_region)

#Back void gap
back_void_gap_region = -outer_cylinder & +plug_back & -tally_cell_front
back_void_gap_cell = openmc.Cell(region=back_void_gap_region)

#Tally 1
tally_1_region = -tally_cell_1_cylinder & +tally_cell_front & -boundary_back
tally_1_cell = openmc.Cell(region=tally_1_region)

#Tally 2
tally_2_region = -tally_cell_2_cylinder & +tally_cell_1_cylinder & +tally_cell_front & -boundary_back
tally_2_cell = openmc.Cell(region=tally_2_region)

#Tally 3
tally_3_region = -tally_cell_3_cylinder & +tally_cell_2_cylinder & +tally_cell_front & -boundary_back
tally_3_cell = openmc.Cell(region=tally_3_region)

#Tally 4
tally_4_region = -tally_cell_4_cylinder & +tally_cell_3_cylinder & +tally_cell_front & -boundary_back
tally_4_cell = openmc.Cell(region=tally_4_region)

#Tally cells
tally_cells_region = -outer_cylinder & +tally_cell_4_cylinder & +tally_cell_front & -boundary_back
tally_cells_cell = openmc.Cell(region=tally_cells_region)

#Boundary cell
boundary_region = -boundary_front & +boundary_back & +outer_cylinder
boundary_cell = openmc.Cell(region=boundary_region)

universe = openmc.Universe(0,name='Main Universe',cells=[plug_hole_cell,SW_cylinder_cell,steel_back_plate_cell,inner_cylinder_cell,gap_2cm_cell,steel_cylinder_cell,source_cell,front_void_gap_cell,back_void_gap_cell,tally_1_cell,tally_2_cell,tally_3_cell,tally_4_cell,tally_cells_cell,boundary_cell])

geom = openmc.Geometry(universe)
geom.export_to_xml()
print ("geo export complete")

############
# Settings #
############
#Initialise openmc settings
settings = openmc.Settings()
#Initialise mesh limits and dimensions
lower_left = (-100,110,-100)
upper_right = (100,660,100)
pitch = 5
wkDir = '/home/mmavis/openmc/openmc3/bin'
volCalc = openmc.calculate_voxel_volumes()
volCalc.lower_left = lower_left
volCalc.upper_right = upper_right
volCalc.pitch = pitch
bounds = volCalc.calculateBounds()
bounds_x = bounds.xbounds
bounds_y = bounds.ybounds
bounds_z = bounds.zbounds
nint = bounds.nint
print(nint)

print(wkDir + "/openmc")

print("defining source parameters")


# Makes the 3d "cube" style geometry 
vox_plot = openmc.Plot()
vox_plot.type = 'voxel'
vox_plot.width = (300., 2000., 300.)
vox_plot.pixels = (300, 300, 300)
vox_plot.filename = 'plot_3d_vol_test'
vox_plot.color_by = 'material'
plots = openmc.Plots([vox_plot])
plots.export_to_xml()
print("plot exported")


#Run OpenMC
print("running OpenMC")
openmc.plot_geometry(openmc_exec= wkDir + "/openmc")


#os.system('/home/mmavis/openmc/scripts/openmc-voxel-to-silovtk plot_3d_vol_test.h5 -o Iter.vti')


##########
# Source #
##########
# Define source parameters
batches = 100
#settings.verbosity = 9
settings.output = {'tallies': False}
settings.batches = batches
settings.inactive = 9
particles = 1e8
settings.particles = int(particles)
print("Particle = " + str(particles))
settings.run_mode = 'fixed source'
histories = (settings.particles*batches) - (settings.inactive*batches)
source = openmc.Source()
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])
#source.energy = openmc.stats.Muir(e0=14080000.0, m_rat=5.0, kt=20000.0)
strength = 2e19
source.strength = int(strength)
source.space = openmc.stats.Box((-100,0,-100),(100,10,100))
settings.source = source
xml_source = source.to_xml_element()
settings.export_to_xml()
print("settings exported")
print("defining tally mesh")

#Define Tally Mesh parameters
mesh_upper_right = upper_right
mesh_lower_left = lower_left
mesh = openmc.RegularMesh()
xdimension = bounds.nint[0] -1
ydimension = bounds.nint[1] -1
zdimension = bounds.nint[2] -1
mesh.dimension = [xdimension, ydimension, zdimension]
mesh.width = [pitch,pitch,pitch]
mesh.lower_left = mesh_lower_left


#Tallies
tallies = openmc.Tallies()

#Mesh Filter
mesh_filter = openmc.MeshFilter(mesh)
mesh_tally = openmc.Tally(1,name='tallies_on_mesh')

#Energy Filter
energy_bins = openmc.mgxs.GROUP_STRUCTURES['VITAMIN-J-175']
energy_filter = openmc.EnergyFilter(energy_bins)

mesh_tally.filters = [mesh_filter] #energy_filter
tallies_to_plot = 'flux'
mesh_tally.scores = [tallies_to_plot]
tallies.append(mesh_tally)
tallies.export_to_xml()
print("Tallies exported")
print ("defining plot")



#Run OpenMC
print("running OpenMC")
openmc.run(openmc_exec= wkDir + '/openmc',threads=16,mpi_args=['mpiexec', '-n', '12'])

sp = openmc.StatePoint('statepoint.'+str(batches)+'.h5')
tally = sp.get_tally(scores=['flux'])
print(tally.mean.shape)
print(tally)
flux = tally.get_slice(scores=['flux'])
print(flux.mean.shape)
flux = flux/125
flux.mean.shape = (xdimension, ydimension, zdimension)

flux_mean_data_1 = (flux.mean[20,:,:])
plt.imshow(flux_mean_data_1,cmap=cm.jet, norm=LogNorm(vmin=1e8, vmax=1e15))
plt.colorbar()
plt.savefig('X.png')
plt.clf()

flux_mean_data_2 = (flux.mean[:,55,:])
plt.imshow(flux_mean_data_2,cmap=cm.jet, norm=LogNorm(vmin=1e8, vmax=1e15))
plt.colorbar()
plt.savefig('Y.png')
plt.clf()

flux_mean_data_3 = (flux.mean[:,:,20])
plt.imshow(flux_mean_data_3,cmap=cm.jet, norm=LogNorm(vmin=1e8, vmax=1e15))
plt.colorbar()
plt.savefig('Z.png')
plt.clf()

#print ("defining volumes to calculate")
volumes = volCalc.calculateVolumes(wkDir, mats)


#Following is the code for taking the results read in from the volume_*.h5 files and outputting the results to the posmat.txt file.

####################
# File Output Code #
####################
posmat = openmc.posmat()
posmat.generator(volumes,bounds)

meshtal = openmc.meshtal()
meshtal.generator(wkDir,bounds,settings,mesh_tally,energy_bins)