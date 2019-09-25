#!/usr/bin/env python3

import openmc
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np
import datetime
from string import Formatter
from sys import getsizeof

class meshtal():
    def generator(self,wkDir,bounds,settings,mesh_tally,energy_bins):
        #Create meshtal file
        histories = (settings.particles*settings.batches) - (settings.inactive*settings.batches)
        sp = openmc.StatePoint(wkDir + '/statepoint.'+str(settings.batches)+'.h5')
        flux_tally = sp.get_tally(scores=mesh_tally.scores)
        #Create a dataframe from the tally results
        tally_results = flux_tally.get_pandas_dataframe()
        tally_number = 1 #figure out a way to change this to array of ids
        #Remove unnecessary data
        tally_results.drop("energy high [eV]",axis=1, inplace=True)
        tally_results.drop("nuclide",axis=1, inplace=True)
        tally_results.drop("score",axis=1, inplace=True)
        df_length = len(tally_results.index)

        #Turn the mesh index into actual xyz coordinates of the midpoints of the mesh cells
        n=1
        for n in range(bounds.nint[0]):
            tally_results[('mesh 1','x')].replace(n,(bounds.lower_left[0]+(bounds.pitch*(n-1)+(bounds.pitch/2))),inplace=True)

        n=1
        for n in range(bounds.nint[1]):
            tally_results[('mesh 1','y')].replace(n,(bounds.lower_left[1]+(bounds.pitch*(n-1)+(bounds.pitch/2))),inplace=True)                

        n=1
        for n in range(bounds.nint[2]):
            tally_results[('mesh 1','z')].replace(n,(bounds.lower_left[2]+(bounds.pitch*(n-1)+(bounds.pitch/2))),inplace=True)

        #Resructure the data frame into the right format
        tally_results.columns = [' '.join(col).strip() for col in tally_results.columns.values]
        tally_results.rename(columns={'mesh 1 x': 'X','mesh 1 y': 'Y','mesh 1 z': 'Z'},inplace=True)

        tally_results['mean'] = tally_results['mean']/1e6
        tally_results['energy low [eV]'] = tally_results['energy low [eV]']/1e6
        if 'std. dev.' in tally_results.columns:
            tally_results = tally_results[['energy low [eV]','X','Y','Z','mean','std. dev.']]
            tally_results.head()
            tally_results['std. dev.'] = tally_results['std. dev.']/1e6
        else:
            tally_results = tally_results[['energy low [eV]','X','Y','Z','mean']]
            tally_results.head()
            n=0
            std_dev_zeros = []
            for n in range(df_length):
                std_dev_zeros.append(0)
            tally_results['std. dev.'] = std_dev_zeros
        #Rename columns
        tally_results.rename(columns={'energy low [eV]': 'Energy', 'mean':'Results', 'std. dev.': 'Rel Error'}, inplace=True)
        #Sort the values by energy x,y,z
        tally_results.sort_values(by=['Energy','X','Y','Z'],ascending=[True,True,True,True],inplace=True)
        #Remove index
        tally_results.reset_index(inplace=True,drop=True)

        energy_bins_np = np.array([energy_bins])
        energy_bins_np = energy_bins_np/1e6
        #Turn dataframe into numpy array
        tally_results = tally_results.to_numpy()
        #Turn titles for columns into an array
        tally_titles = np.array([['Energy','X','Y','Z','Results','Rel Error']])

        #Open the posmat.txt file and clear it if it already exists. If not create it.
        open('meshtal.msht', 'w').close()
        #Open the meshtal file in append mode
        with open('meshtal.msht','a+') as tally_file:
            #Set todays date
            today = datetime.datetime.now()
            #Write file header
            tally_file.write("OpenMC   version 0.11.0  "+today.strftime('%d/%m/%Y')+ " probid =  "+today.strftime('%d/%m/%Y %H:%M:%S'+"\n"))
            #Write mesh type
            tally_file.write(" Cube"+"\n")
            #Write histories
            tally_file.write(" Number of histories used for normalizing tallies ="+"{0:17.2f}".format(histories)+"\n"+"\n")
            #Write mesh tally number
            tally_file.write(" Mesh Tally Number"+"{0:10d}".format(tally_number)+"\n"+" neutron  mesh tally." +"\n"+"\n")
            #Write mesh tally boundries
            tally_file.write(" Tally bin boundaries:"+"\n")
            #Output X bounds
            tally_file.write("    X direction:")
            for item in bounds.xbounds:
                tally_file.write('{0:10.2f}'.format(item))
            tally_file.write("\n")
            #Output Y bounds
            tally_file.write("    Y direction:")
            for item in bounds.ybounds:
                tally_file.write('{0:10.2f}'.format(item))
            tally_file.write("\n")
            #Output Z bounds
            tally_file.write("    Z direction:")
            for item in bounds.zbounds:
                tally_file.write('{0:10.2f}'.format(item))
            tally_file.write("\n")
            #Write energy bounds
            tally_file.write("    Energy bin boundaries: ")
            np.savetxt(tally_file,(energy_bins_np),fmt="%8.2e",newline='\n')
            tally_file.write("\n")
            #Write column titles
            np.savetxt(tally_file,tally_titles,fmt="%11s %9s %9s %9s %11s %11s ",delimiter=',',newline='\n')
            #Write tally data
            np.savetxt(tally_file,tally_results,fmt="%11.3e %9.3f %9.3f %9.3f %11.5e %11.5e",delimiter=' ',newline='\n')
        return