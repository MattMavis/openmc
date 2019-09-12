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
        lower_left_x_int, lower_left_y_int, lower_left_z_int = self.lower_left
        upper_right_x_limit, upper_right_y_limit, upper_right_z_limit = self.upper_right
        pitch = self.pitch
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

        self.xbounds = bounds_x
        self.ybounds = bounds_y
        self.zbounds = bounds_z
        nint = (len(bounds_x),len(bounds_y),len(bounds_z))
        self.nint = nint
        return (self)


    def calculateVolumes(self, wkDir, mats):
        settings = openmc.Settings()
        #Initialise Coording array to store mid co-ordinates of the cells
        #output bounds
        # Work out how many cells in each direction. These are the dimensions for the meshtal file
        pitch = self.pitch
        bounds = self.calculateBounds()
        resultCoord = []
        results = []
        j = 0
        b=0
        #Loop until all mesh has been 'scan' and volume data is produced for all cells in the mesh
        lower_left_x = self.lower_left[0]

        while lower_left_x<= self.upper_right[0]:
            lower_left_y = self.lower_left[1]
            b +=1
            print("Starting Batch " + str(b))
            #Initialise vol_calc array to store volume data for each cell of the mesh
            vol_calc_array = []

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

                    upper_right_x = lower_left_x + pitch
                    upper_right_y = lower_left_y + pitch
                    upper_right_z = lower_left_z + pitch
                    upper_right = (upper_right_x, upper_right_y, upper_right_z)
                    mid_coord = (lower_left_x + (pitch/2),lower_left_y + (pitch/2),lower_left_z + (pitch/2))
                    resultCoord.append(mid_coord)
                    vol_calc = openmc.VolumeCalculation(mats, 10000, lower_left, upper_right)
                    vol_calc_array.append(vol_calc)
                    j += 1
                    lower_left_z += pitch

                lower_left_y += pitch

            lower_left_x += pitch

            settings.volume_calculations = vol_calc_array
            settings.export_to_xml()
            #Run Openmc in volume calculation mode
            print("calculating volumes")
            openmc.calculate_volumes(openmc_exec= wkDir + '/openmc',mpi_args=['mpiexec', '-n', '16'],threads=16)
            print("Reading Results ...")
            m=1
            while m <= ((len(bounds.ybounds)-1)*(len(bounds.zbounds)-1)):
                vol_calc.load_results(wkDir + '/volume_' + str(m) +'.h5')
                results.append(vol_calc.volumes)
                m+=1
            print(str(m-1) + " Results Read!")
            print("Cleaning Working Directory ...")
            os.system('rm ' + wkDir + '/volume_*.h5')
            print("Directory Cleaned!")
            vol_calc_array.clear()
            self.results = results
            self.resultCoord = resultCoord
            self.numResults = j
        
        return (self)