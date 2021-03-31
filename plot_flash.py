#This code is for making plots from FLASH output files
#Created by Malia Jenks

import h5py
import numpy as np
import pylab
from math import sqrt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

import os 
import re 
import sys 
import optparse

def read_option():
    usage = "usage: [%prog] [options]  Flash_outputs \n" 
    usage += "Takes all the outputs and makes plots." 

    parser = optparse.OptionParser(usage=usage)

    parser.add_option("-d", "--output_path",
                      dest = "output_path",
                      type = "string",
                      help =  "Directory for the output files. Default is plots.",
                      default = "plots/")


    option,args = parser.parse_args()
    option.args = args
    
    if not option.output_path.endswith("/"):
        option.output_path += "/"

    if not os.path.exists(option.output_path):
        os.makedirs(option.output_path)

    return option 


if __name__ ==  "__main__": 
    option = read_option() 

    for flash_output in option.args:

        filename = flash_output
        output_path = option.output_path
        index = flash_output[-4:]


        #Read in hdf5 file
        h5file = h5py.File(filename, 'r')

        #get the number of blocks and the dimensions of each block
        numblocks = np.shape(h5file['temp'])[0]
        blk_r = np.shape(h5file['temp'])[3]
        blk_z = np.shape(h5file['temp'])[2]

        #First find the r and z values of the centers of blocks
        #And the r and z ranges of the total grid
        x = []
        y = []
        rmin, rmax, = h5file['bounding box'][0][0] #prefilled with the values from block 0
        zmin, zmax, = h5file['bounding box'][0][1]
        for i in range(0,numblocks):
            x.append(h5file['coordinates'][i][0])
            y.append(h5file['coordinates'][i][1])
            if h5file['bounding box'][i][0][0] < rmin:
                rmin = h5file['bounding box'][i][0][0]
            if h5file['bounding box'][i][0][1] > rmax:
                rmax = h5file['bounding box'][i][0][1]
            if h5file['bounding box'][i][1][0] < zmin:
                zmin = h5file['bounding box'][i][1][0]
            if h5file['bounding box'][i][1][1] > zmax:
                zmax = h5file['bounding box'][i][1][1]
            
        r_centers = np.sort(np.unique(x))
        z_centers = np.sort(np.unique(y))
        r_centers = r_centers.tolist()
        z_centers = z_centers.tolist()
        
        #Make the r and z arrays
        rstep = (rmax - rmin)/(sqrt(numblocks)*blk_r)
        zstep = (zmax - zmin)/(sqrt(numblocks)*blk_z)

        r = []
        z = []

        #I am making the r and z arrays together because they are the same size
        #if this changes this will need to be split to two loops
        for j in range(int(sqrt(numblocks)*blk_r)):
            r.append(rmin + rstep*j)
            z.append(zmin + zstep*j)

        r = np.asarray(r)
        z = np.asarray(z)

        #Make empty numpy arrays for coordinates and variables
        var_shape = (len(r),len(z))
        temp = np.zeros(var_shape)
        dens = np.zeros(var_shape)
        vz = np.zeros(var_shape)
        vr = np.zeros(var_shape)
        vph = np.zeros(var_shape)

        #Fill in the variable arrays

        for n in range(numblocks):
            center = h5file['coordinates'][n]
            r_ind = r_centers.index(center[0])
            z_ind = z_centers.index(center[1])
            r_start = r_ind*blk_r
            z_start = z_ind*blk_z
            temp0 = np.asarray(h5file['temp'][n][0])
            dens0 = np.asarray(h5file['dens'][n][0])
            vz0 = np.asarray(h5file['vely'][n][0])
            vr0 = np.asarray(h5file['velx'][n][0])
            vph0 = np.asarray(h5file['velz'][n][0])
            for l in range(blk_r):
                r0 = r_start + l
                for m in range(blk_z):
                    z0 = z_start + m
                    temp[r0,z0] = temp0[m,l]
                    dens[r0,z0] = dens0[m,l]
                    vz[r0,z0] = vz0[m,l]
                    vr[r0,z0] = vr0[m,l]
                    vph[r0,z0] = vph0[m,l]

        h5file.close()
                    
        #The variable arrays got transposed by hdf5 so correct that
        temp = np.transpose(temp)
        dens = np.transpose(dens)
        vz = np.transpose(vz)
        vr = np.transpose(vr)
        vph = np.transpose(vph)
        
        vphi = np.fabs(vph)


        #make plots

        golden = (1.0 + sqrt(5.0)) / 2.0
 
        figprops = dict(figsize=(8., 8./golden), dpi=128)
        adjustprops = dict(left=0.1, bottom=0.12, right=0.97, top=0.93, wspace=0.3, hspace=0.3)
 
        fig1 = pylab.figure(**figprops)
        fig1.subplots_adjust(**adjustprops)
 
        ax = fig1.add_subplot(111)
 
        cax = ax.pcolormesh(r, z, temp, norm=LogNorm(vmin=5e7, vmax=5e9), cmap='OrRd')
 
        ax.set_xlabel("r", fontsize = 20)
        ax.set_ylabel("z", fontsize = 20)
        ax.set_xlim([0,8e9])
        ax.set_ylim([-4e9,4e9])
        fig1.colorbar(cax)
        fig1.suptitle('Temperature', fontsize = 20)
        name = output_path + 'cyl_grid_temp' + index
        pngname = name + ".png"
        fig1.savefig(pngname)
 
        pylab.close
 
        fig2 = pylab.figure(**figprops)
        fig2.subplots_adjust(**adjustprops)
 
        ax1 = fig2.add_subplot(111)
 
        cax = ax1.pcolormesh(r, z, dens, norm=LogNorm(vmin=1e2, vmax=1e9), cmap='PuRd')
 
        ax1.set_xlabel("r", fontsize = 20)
        ax1.set_ylabel("z", fontsize = 20)
        ax1.set_xlim([0,8e9])
        ax1.set_ylim([-4e9,4e9])
        fig2.colorbar(cax)
        fig2.suptitle('Density', fontsize = 20)
        name = output_path + 'cyl_grid_dens' + index
        pngname = name + ".png"
        fig2.savefig(pngname)
 
        pylab.close
 
        fig3 = pylab.figure(**figprops)
        fig3.subplots_adjust(**adjustprops)
 
        ax1 = fig3.add_subplot(111)
 
        cax = ax1.pcolormesh(r, z, vphi, norm=Normalize(vmin=0, vmax=5e8), cmap='BuPu')
#        cax = ax1.pcolormesh(r, z, vphi, norm=LogNorm(vmin=np.amin(vphi), vmax=np.amax(vphi)), cmap='BuPu')
 
        ax1.set_xlabel("r", fontsize = 20)
        ax1.set_ylabel("z", fontsize = 20)
        fig3.colorbar(cax)
        fig3.suptitle('Absolute Value of Phi Velocity', fontsize = 20)
        name = output_path + 'cyl_grid_vphi' + index
        pngname = name + ".png"
        fig3.savefig(pngname)
 
        pylab.close
 
        fig4 = pylab.figure(**figprops)
        fig4.subplots_adjust(**adjustprops)
 
        ax1 = fig4.add_subplot(111)
 
        cax = ax1.pcolormesh(r, z, vz, norm=Normalize(vmin=-5e8,vmax=5e8), cmap='Spectral')
 
        ax1.set_xlabel("r", fontsize = 20)
        ax1.set_ylabel("z", fontsize = 20)
        fig4.colorbar(cax)
        fig4.suptitle('Vz', fontsize = 20)
        name = output_path + 'cyl_grid_vz' + index
        pngname = name + ".png"
        fig4.savefig(pngname)
 
        pylab.close
 
        fig5 = pylab.figure(**figprops)
        fig5.subplots_adjust(**adjustprops)
 
        ax1 = fig5.add_subplot(111)
 
        cax = ax1.pcolormesh(r, z, vr, norm=Normalize(vmin=0,vmax=5e8), cmap='copper')
 
        ax1.set_xlabel("r", fontsize = 20)
        ax1.set_ylabel("z", fontsize = 20)
        fig5.colorbar(cax)
        fig5.suptitle('Vr', fontsize = 20)
        name = output_path + 'cyl_grid_vr' + index
        pngname = name + ".png"
        fig5.savefig(pngname)
 
        pylab.close

        pylab.close('all')



