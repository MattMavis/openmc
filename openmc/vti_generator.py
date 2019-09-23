#!/usr/bin/env python3

import openmc
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np
import pandas as pd
import datetime
import vtk
import linecache
from string import Formatter
import sys

class vti_generator(object):

    def progressBar(value, endvalue, bar_length=20):

        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    def generator(self, wkDir,batches):
        # Get the total amount of particles from the MCR2S commonfile
        total_source = linecache.getline(wkDir+'/commonfile',5)
        total_source_str, Total_particles = total_source.split(" ")
        Total_particles = float(Total_particles)
        # Open the statepoint file
        filename = wkDir + '/statepoint.'+str(batches)+'.h5'
        f = h5py.File(filename,'r')
        # Get the number of realisations (batches)
        n_realisations = f['n_realizations'][()]
        # Get the results from the statepoint file from the tally
        Tally_1 = f['/tallies/tally 1']
        results = f['/tallies/tally 1/results'][()]
        print(results)
        # Seperate the results into the flux and the std dev
        result_1, result_2 = results.T
        # Reshapes them into columns
        result_1 = np.reshape(result_1,(-1,1))
        result_2 =np.reshape(result_2,(-1,1))
        #Calculate the mean flux by dividing by the realisations and multiplying by the total particles
        mean = (result_1/n_realisations)*Total_particles
        #Do the same to the std dev
        l=0
        std_dev = []
        print(mean.size)
        std_dev = (((result_2/n_realisations)-(result_1)**2)/(n_realisations-1))
        std_dev = ((np.sqrt(std_dev))*Total_particles)
        # Get the mesh dimensions
        mesh_dimensions = f['/tallies/meshes/mesh 1/dimension'][()]
        # Seperate the mesh dimension into the x,y,z
        nx, ny, nz = mesh_dimensions.T
        # Get the mesh width
        mesh_width = f['/tallies/meshes/mesh 1/width'][()]
        # Seperate into x,y,z
        x_width, y_width, z_width = mesh_width.T
        # Create an array with the widths
        width = (x_width,y_width,z_width)
        # Get the lower left co-ordinates of the mesh
        mesh_lower_left = f['/tallies/meshes/mesh 1/lower_left'][()]
        print(mesh_lower_left)
        # Seperate the co-ordinates
        x_lower_left, y_lower_left, z_lower_left = mesh_lower_left.T
        lower_left = (x_lower_left, y_lower_left, z_lower_left)

        # Redundent?
        # Produce a dataframe with the mean flux
        energy_df = pd.DataFrame({"Mean": mean[:,0]})
        # Turn the data frame into an array
        array = energy_df.to_numpy()

        # Close the file
        f.close()

        # Initialise vtkMultiBlockDataSet
        blocks = vtk.vtkMultiBlockDataSet()
        blocks.SetNumberOfBlocks(5)
        block_idx = 0

        # Initialise vtkImageData with the mesh dimensions
        grid = vtk.vtkImageData()
        grid.SetDimensions(nx+1, ny+1, nz+1)
        grid.SetOrigin(*lower_left)
        grid.SetSpacing(*width)

        print("Getting mesh mean flux data...")
        # Initialise vtkDoubleArray so it can recieve the mean flux data
        data = vtk.vtkDoubleArray()
        data.SetName("Flux")
        s = 0
        pos = []
        cell_flux = []
        # Loop over the mesh cell in the order z,y,x and assign there flux value
        for x in range(nx):
            vti_generator.progressBar(x,nx)
            for y in range(ny):
                for z in range(nz):
                    cell_flux = np.append(cell_flux, mean[s])
                    i = z*nx*ny + y*nx + x
                    pos = np.append(pos,i)
                    s+=1          
        value = 0
        # Loop over the flux data and set it to be in the order needed to add it to the cells
        for p in range(nx*ny*nz):
            data.InsertNextValue(cell_flux[p])
            value += 1
        print(value)
        # Set the grid cell data
        grid.GetCellData().AddArray(data)
        print(data)
        del data

        blocks.SetBlock(block_idx, grid)
        block_idx += 1
        print("Getting mesh std_dev flux data ...")
        # Do the same with standard deviation
        StdDevData = vtk.vtkDoubleArray()
        StdDevData.SetName("StdDev")
        s = 0
        pos = []
        StdDeV = []
        for x in range(nx):
            vti_generator.progressBar(x,nx)
            for y in range(ny):
                for z in range(nz):
                    StdDeV = np.append(StdDeV, std_dev[s])
                    i = z*nx*ny + y*nx + x
                    pos = np.append(pos,i)
                    s+=1
            print(str(x+1) + " out of " + str(nx+1) + " batches processed")
        for p in range(nx*ny*nz):
            StdDevData.InsertNextValue(StdDeV[p])
            value += 1
        print(value)
        grid.GetCellData().AddArray(StdDevData)
        print(StdDevData)
        del StdDevData
        blocks.SetBlock(block_idx, grid)
        block_idx += 1

        # Set writer to write the multiblock data
        writer = vtk.vtkXMLMultiBlockDataWriter()
        # Check version of VTK and use the right function for the version
        if vtk.vtkVersion.GetVTKMajorVersion() > 5:
            writer.SetInputData(blocks)
        else:
            writer.SetInput(blocks)
        writer.SetFileName("Flux.vti")
        # Writes the vtk/vti file
        writer.Write()
        return







