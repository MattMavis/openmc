#!/usr/bin/env python3
import h5py
import numpy as np
import pandas as pd
import openmc
from pyevtk.hl import gridToVTK 
import random as rand

lower_left_x_int = 0
lower_left_y_int = -523
lower_left_z_int = -523
upper_right_x_limit = 523 
upper_right_y_limit = 0
upper_right_z_limit = 523
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

filename = 'statepoint.2.h5'
f = h5py.File(filename,'r')
print("Keys: %s" % f.keys())
n_realisations = f['n_realizations'][()]
print(n_realisations)
Tally_1 = f['/tallies/tally 1']
print(Tally_1)
print("Keys: %s" % Tally_1.keys())
results = f['/tallies/tally 1/results'][()]
print(results)
batch_1, batch_2 = results.T
batch_1 = np.reshape(batch_1,(-1,1))
batch_2 =np.reshape(batch_2,(-1,1))

energy_bins = f['/tallies/filters/filter 2/n_bins'][()]
print(energy_bins)

x = []
y = []
z = []
mesh_dimensions = f['/tallies/meshes/mesh 1/dimension'][()]
print(mesh_dimensions)
x_dimension, y_dimension, z_dimension = mesh_dimensions.T
print(x_dimension)
print(y_dimension)
print(z_dimension)
mesh_width = f['/tallies/meshes/mesh 1/width'][()]
print(mesh_width)
x_width, y_width, z_width = mesh_width.T
print(x_width)
print(y_width)
print(z_width)
mesh_lower_left = f['/tallies/meshes/mesh 1/lower_left'][()]
print(mesh_lower_left)
x_lower_left, y_lower_left, z_lower_left = mesh_lower_left.T
print(x_lower_left)
print(y_lower_left)
print(z_lower_left)

for i in range(x_dimension):
    for j in range(energy_bins):
        x = np.append(x,x_lower_left + i*x_width)
print(x)
for i in range(y_dimension):
    for j in range(energy_bins):
        y = np.append(y,y_lower_left + i*y_width)
print(y)
for i in range(z_dimension):
    for j in range(energy_bins):
        z = np.append(z,z_lower_left + i*z_width)
print(z)


energy_df = pd.DataFrame({"Mean": batch_1[:,0]/n_realisations,"Std. Dev": batch_2[:,0]/n_realisations})
print(energy_df)
print(energy_df.loc[[12]])
batches = 2
tallies_to_plot = 'flux'




""" sp = openmc.StatePoint('statepoint.'+str(batches)+'.h5')
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
print(tally_results) """
