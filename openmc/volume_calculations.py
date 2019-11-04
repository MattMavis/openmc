from collections import OrderedDict
from collections.abc import Iterable, Mapping
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import warnings
import os
import numpy as np
import pandas as pd
import h5py
from uncertainties import ufloat

import openmc
import openmc.checkvalue as cv

_VERSION_VOLUME = 1

class calculate_voxel_volumes(object):

    def calculateBounds(self):
        #Set the input bounds of the mesh. These bounds are the lower left and upper right corner. The width 
        #of the mesh cell and the starting mesh cell position.
        lower_left_x_int, lower_left_y_int, lower_left_z_int = self.lower_left
        upper_right_x_limit, upper_right_y_limit, upper_right_z_limit = self.upper_right
        width = self.width
        xbound_value = lower_left_x_int
        ybound_value = lower_left_y_int
        zbound_value = lower_left_z_int

        bounds_x = [xbound_value]
        bounds_y = [ybound_value]
        bounds_z = [zbound_value]
        
        #Iterate over the value of the x,y,z bounds until it exceeds the upper right x,y,z limits and store those value
        #in the bounds arrays.
        while xbound_value < upper_right_x_limit:
            xbound_value += width
            bounds_x.append(xbound_value)

        while ybound_value < upper_right_y_limit:
            ybound_value += width
            bounds_y.append(ybound_value)

        while zbound_value < upper_right_z_limit:
            zbound_value += width
            bounds_z.append(zbound_value)
        
        #Append bounds arrays to self, along with the amount of bounds in each direction (x,y,z), so that it can be returned
        #to Main
        self.xbounds = bounds_x
        self.ybounds = bounds_y
        self.zbounds = bounds_z
        nint = (len(bounds_x),len(bounds_y),len(bounds_z))
        self.nint = nint
        self.width = width
        return (self)


    def calculateVolumes(self, wkDir, mats):
        
        settings = openmc.Settings()
        width = self.width
        #Call the function calculateBounds to calcualted the bounds of the mesh.
        bounds = self.calculateBounds()
        
        #If the amount of mesh cell exceeds 15000, run in a batch mode to reduce memory and disk space overhead needed. 
        #Else run as normal
        if (((bounds.nint[0]-1)*(bounds.nint[1]-1)*(bounds.nint[2]-1)) >= 15000):
            #Initialise arrays
            #resultCoord stores the mid point co-ordinate of each mesh cell
            #results stores the volume data generated from the stochastic volume calculation for each mesh cell
            #vol_calc_array contains the settings for each volume calculation for each mesh cell
            resultCoord = []
            results = []
            vol_calc_array = []
            #Set loop and batch number to 0 (Both varaiable are used for debugging)
            j = 0
            b = 0
            #Loop until all mesh has been 'scan' and volume data is produced for all cells in the mesh.
            
            #Reset lower left x co-ordinate to its initial position
            lower_left_x = self.lower_left[0]
            while lower_left_x<= self.upper_right[0]:
                #Reset lower left y co-ordinate to its initial position
                lower_left_y = self.lower_left[1]
                b += 1
                print("Starting Batch " + str(b))
                
                #IF lower left x coordinate exceeds or equals the upper right x limit coordinate terminate the loop
                if lower_left_x >= self.upper_right[0]:
                            break

                while lower_left_y<= self.upper_right[1]:
                    #Reset lower left z co-ordinate to its initial position
                    lower_left_z = self.lower_left[2]
                    
                    #IF lower left y coordinate exceeds or equals the upper right y limit coordinate terminate the loop
                    if lower_left_y >= self.upper_right[1]:
                            break            

                    while lower_left_z<= self.upper_right[2]:
                        #Set the lower left coordinate of the bounding box of interest to the values for x,y,z.
                        lower_left = (lower_left_x, lower_left_y, lower_left_z)
                        
                        #IF lower left z coordinate exceeds or equals the upper right z limit coordinate terminate the loop
                        if lower_left_z >= self.upper_right[2]:
                            break
                        #Work out the upper bounds of the bounding box of interest
                        upper_right_x = lower_left_x + width
                        upper_right_y = lower_left_y + width
                        upper_right_z = lower_left_z + width
                        upper_right = (upper_right_x, upper_right_y, upper_right_z)
                        #Find mid co-ordinate of selected voxel and add them to the resultCoord array.
                        mid_coord = (lower_left_x + (width/2),lower_left_y + (width/2),lower_left_z + (width/2))
                        resultCoord.append(mid_coord)
                        #Create settings for a volume calculation in selected mesh by supplying the domains to look for,e.g materials,
                        #the amount of samples, the lower left coordinate and upper right coordinate of the bounding box.
                        vol_calc = openmc.VolumeCalculation(mats, 10000, lower_left, upper_right)
                        #Add volume calculation setting to array of volume calculations
                        vol_calc_array.append(vol_calc)
                        j += 1
                        lower_left_z += width

                    lower_left_y += width

                lower_left_x += width
                
                #Save the volume calculation settings for all mesh cells to openMC's settings variable and export them to the XML
                #file
                settings.volume_calculations = vol_calc_array
                settings.export_to_xml()
                #Run Openmc in volume calculation mode
                print("calculating volumes")
                openmc.calculate_volumes(openmc_exec= wkDir + '/openmc',mpi_args=['mpiexec', '-n', '16'],threads=16)
                print("Reading Results ...")
                m=1
                print(j)
                while m <= ((len(bounds.ybounds)-1)*(len(bounds.zbounds)-1)):
                    #Read in volume results from volume files
                    vol_calc.load_results(wkDir + '/volume_' + str(m) +'.h5')
                    #Add volume results to results array
                    results.append(vol_calc.volumes)
                    m+=1
                print(str(m-1) + " Results Read!")
                #Clean the working directory by deleting work files and free the memory the arrays were using.
                print("Cleaning Working Directory ...")
                os.system('rm ' + wkDir + '/volume_*.h5')
                print("Directory Cleaned!")
                vol_calc_array.clear()
            self.results = results
            self.resultCoord = resultCoord
            self.numResults = j
            results.clear()
            resultCoord.clear()
            
        else:
            resultCoord = []
            results = []
            vol_calc_array = []
            j = 0
            b=0
            #Loop until all mesh has been 'scan' and volume data is produced for all cells in the mesh
            lower_left_x = self.lower_left[0]

            if lower_left_x >= self.upper_right[0]:
                            break

                while lower_left_y<= self.upper_right[1]:
                    lower_left_z = self.lower_left[2]

                    if lower_left_y >= self.upper_right[1]:
                            break            

                    while lower_left_z<= self.upper_right[2]:
                        lower_left = (lower_left_x, lower_left_y, lower_left_z)
                        
                        if lower_left_z >= self.upper_right[2]:
                            break

                        upper_right_x = lower_left_x + width
                        upper_right_y = lower_left_y + width
                        upper_right_z = lower_left_z + width
                        upper_right = (upper_right_x, upper_right_y, upper_right_z)
                        mid_coord = (lower_left_x + (width/2),lower_left_y + (width/2),lower_left_z + (width/2))
                        resultCoord.append(mid_coord)
                        vol_calc = openmc.VolumeCalculation(mats, 10000, lower_left, upper_right)
                        vol_calc_array.append(vol_calc)
                        j += 1
                        lower_left_z += width

                    lower_left_y += width

                lower_left_x += width

            settings.volume_calculations = vol_calc_array
            settings.export_to_xml()
            #Run Openmc in volume calculation mode
            print("calculating volumes")
            openmc.calculate_volumes(openmc_exec= wkDir + '/openmc',mpi_args=['mpiexec', '-n', '16'],threads=16)
            print("Reading Results ...")
            print(j)
            m=1
            while (m <= j):
                vol_calc.load_results(wkDir + '/volume_' + str(m) +'.h5')
                results.append(vol_calc.volumes)
                m+=1
            print(str(m-1) + " Results Read!")
            #Clean the working directory by deleting work files and free the memory the arrays were using.
            print("Cleaning Working Directory ...")
            os.system('rm ' + wkDir + '/volume_*.h5')
            print("Directory Cleaned!")
            self.results = results
            self.resultCoord = resultCoord
            self.numResults = j
            vol_calc_array.clear()
            results.clear()
            resultCoord.clear()
        
            return (self)
