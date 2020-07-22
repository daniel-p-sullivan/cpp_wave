"""Plotter built for CFD_Cpp"""
import os
import subprocess
from glob import glob
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

subprocess.call(['make'])			#make the cpp program
subprocess.call(['CFD_CPP_2D.exe'])	#run the cpp program

plt.ioff()

FILES = glob("./out/****.txt")		#glob all the .txts
FILES.sort()						#not really necessary but w/e

FIG = plt.figure()					#initialize the figure
for F in FILES:						#file loop
    AX = FIG.gca(projection='3d')					#mplot3d
    F = open(F, 'r')								#read
    XCELLS = int(F.readline().split("=")[1])		#get XCELLS
    YCELLS = int(F.readline().split("=")[1])		#get YCELLS
    GC = int(F.readline().split("=")[1])			#get GHOSTCELLS
    TIME = float(F.readline().split("=")[1])		#get TIME
    X, Y, H, HU, HV = np.loadtxt(F, unpack=True)	#loadtxt
    F.close()										#close the file asap
    ROWS = 2 * GC + YCELLS				#define ROWS
    COLS = 2 * GC + XCELLS				#define COLS

    X = np.reshape(X, (ROWS, COLS))     #standard data reshape from 1D -> 2D
    Y = np.reshape(Y, (ROWS, COLS))
    H = np.reshape(H, (ROWS, COLS))

    XPLOT = X[GC:-GC, GC:-GC]           #trim off the Ghostcells for plotting
    YPLOT = Y[GC:-GC, GC:-GC]
    HPLOT = H[GC:-GC, GC:-GC]

    AX.plot_surface(XPLOT, YPLOT, HPLOT, rstride=2, cstride=2, linewidth=0.1) #plot
    AX.set_title('Time = %f seconds' % TIME)								  #Time = TIME seconds
    AX.set_zlim(-0.1, 3.5)													  #set the Z limit
    FIG.savefig("%s.png" % os.path.splitext(F.name)[0])						  #save to .png
    print("Wrote %s.png" % os.path.splitext(F.name)[0])						  #console out
    plt.clf()																  #close the fig

os.chdir('./out/')																	 #get into the ./out/ dir
subprocess.call(['ffmpeg', '-i', '%d.png', '-b', '2048k', 'output.avi'])			 #ffmpeg the pngs to avi
subprocess.call(['ffmpeg', '-i', 'output.avi', '-b', '2048k', '-t', '5', 'out.gif']) #ffmpeg the avi to a gif
