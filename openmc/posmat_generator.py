#!/usr/bin/env python3

import openmc
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np
import datetime
from string import Formatter
from sys import getsizeof


class posmat(object):
    def generator(self,volumes,bounds):
        
        k = 0
        l = 0
        outputVolPer = []
        #Calculate the total volume of the sampled areas.
        total_volume = (volumes.pitch**3)                                           
        coord = 'X Cord','Y Cord','Z Cord'                                         
        #Break the results into the 2 seperate arrays.
        mat, vol = zip(*volumes.results[k].items())                                     
        #Open the posmat.txt file and clear it if it already exists. If not create it.
        open('posmat.txt', 'w').close()                                     
        #Open posmat.txt in append mode.
        with open('posmat.txt', 'a+') as file_handler:
            #Output the amount of bounds in each direction
            for item in bounds.nint:
                file_handler.write('{0:>12d}'.format(item))
            file_handler.write("\n")
            #Output X bounds
            file_handler.write("xbounds ")
            for item in bounds.xbounds:
                file_handler.write('{0:10.2f}'.format(item))
            file_handler.write("\n")
            #Output Y bounds
            file_handler.write("ybounds ")
            for item in bounds.ybounds:
                file_handler.write('{0:10.2f}'.format(item))
            file_handler.write("\n")
            #Output Z bounds
            file_handler.write("zbounds ")
            for item in bounds.zbounds:
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
        for k in range(volumes.numResults):
            #Reset l
            l = 0
            #Break the results at array positon k into the 2 seperate arrays.
            mat, vol = zip(*volumes.results[k].items())
            #Find out the amount of values in the results at k.
            numVolumes = len(volumes.results[k])
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
            #Open posmat.txt in append mode.
            with open('posmat.txt', 'a+') as file_handler:
                #Iterate through the values in resultCoord[k] and write them to posmat.txt, each value is in scientific form with 4 decimal places and is justified to the right.
                f=0
                for item in volumes.resultCoord[k]:
                    fmt=["{0:12.4e}","{0:13.4e}","{0:13.4e}"]
                    file_handler.write(fmt[f].format(item))
                    f += 1
                #Iterate through the values in outputVolPer and write them to posmat.txt, each value is a float with 8 places with 6 of them decimal places and is justified to the right.
                for item in outputVolPer:
                    file_handler.write(str(" ")+"{0:8.6f}".format(item))
                #Start new line for next result.
                file_handler.write('\n')
            #Close the file
            file_handler.close()
        return