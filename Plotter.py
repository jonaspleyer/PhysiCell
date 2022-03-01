#!/bin/python3

from matplotlib.animation import FuncAnimation
from matplotlib import cm
import matplotlib.pyplot as plt
import os
import numpy as np
import multiprocessing as mp
import tqdm
from functools import partial

# Import custom modules
from InfoGetter import InfoGetter

class Plotter:
	def __init__(self):
		self.IG = InfoGetter()
		self.style = 'seaborn-deep'
		self.cmap = cm.viridis


	def initializePltSettings(self, sizex=12, sizey=9):
		'''Store the plotting settings for matplotlib and return figure.'''
		plt.style.use(self.style)
		fig = plt.figure(figsize=(sizex, sizey), dpi=100)
		return fig


	def printOutput(self, i, var_name, N_files):
		'''Simple function to show how far in the progress we are.'''
		N = str(len(str(N_files)))
		first = ("[{:0"+ N + "d}").format(i)
		second = ("/{:0" + N + "d}]").format(N_files-1)
		print(first + second + " " + var_name)


	def animate(self, i, substrate_name, ax_sim, fig):
		'''Function used to animate the 2D substrate plots.'''
		M = self.IG.getMicroenv(i)

		self.printOutput(i, substrate_name, len(self.IG.files_mat_micros))

		# Create a grid for plotting later on
		xgrid = M[0, :].reshape(self.IG.numy,self.IG.numx)
		ygrid = M[1, :].reshape(self.IG.numy,self.IG.numx)

		N = self.IG.getSubstrate(i, substrate_name)
		if N.shape[0] != 1:
			raise ValueError("Plotter.animate: Expected 2D microenvironment. Check input files.")

		# Actualy plot the surface
		fig.clear()
		ax_sim = fig.add_subplot(111)


		# Create Titles for plots
		ax_sim.set_title("Substrate " + substrate_name + ": Sim $M(x,y)$")

		# Fill subplots with information
		# ax_sim.plot_surface(xgrid, ygrid, N[0], cmap=self.cmap)
		x = xgrid[0]
		y = ygrid[:,0]
		N_ticks_x = min(len(x), 6 + len(x) % 2)
		N_ticks_y = min(len(y), 6 + len(x) % 2)
		ticks_step_x = int(len(x)/N_ticks_x)
		ticks_step_y = int(len(y)/N_ticks_y)
		ax_sim.set_xticks(ticks_step_x*np.arange(N_ticks_x+1), labels=[x[min(len(x), ticks_step_x*i)] for i in range(N_ticks_x+1)])
		ax_sim.set_yticks(ticks_step_y*np.arange(N_ticks_y+1), labels=[y[min(len(x), ticks_step_y*i)] for i in range(N_ticks_y+1)])
		pos = ax_sim.imshow(N[0], interpolation='bicubic', cmap=self.cmap, vmin=self.N_min, vmax=self.N_max)
		cbar = fig.colorbar(pos, ax=ax_sim, format='%E')#, cax=fig.add_axes([0.1, (mi-self.N_min)/(self.N_max-self.N_min), 0.2, (ma-mi)/(self.N_max-self.N_min)]))
		cbar.minorticks_on()

		# Save individual plots to directory
		fig.savefig("pictures/" + "/plot_substr_" + substrate_name + "_" + "%08d" % i + ".svg")


	def generateSubstrateAnim(self, substrate_name, suffix='gif'):
		'''Creates an animation of the 2D substrate densitites in the simulation steps and stores it (default as gif).'''
		fig = self.initializePltSettings()
		ax_sim = fig.add_subplot(111)#, projection='3d')

		self.__createPictureFolder__("./pictures")

		# Define maximum and minimum value of plotting range for whole plotting time
		N_all = self.IG.getSubstrateSeries(substrate_name)
		self.N_max = max(np.array(N_all).flatten())
		self.N_min = min(np.array(N_all).flatten())

		# Settings for the animation
		fps = 10
		# Animate
		anim = FuncAnimation(fig, self.animate, len(self.IG.files_mat_micros), fargs=[substrate_name, ax_sim, fig], interval=1000/fps)

		# Check if a file already exists and create new if so
		i = 0
		while i < 1000:
			filename = 'substrate_' + substrate_name + "_" + str(i) + '.' + suffix
			if os.path.isfile(filename):
				i += 1
			else:
				print("Start generating animation")
				break

		# Store animation in file
		anim.save(filename, writer='imagemagick')
		print("Saved to " + filename)


	def __createPictureFolder__(self, folder):
		'''Creates a folder if it does not exist.'''
		if not os.path.isdir(folder):
			os.mkdir(folder)


	def exportSinglePicture(self, substrate_name, folder, i):
		fig = self.initializePltSettings()
		ax_sim = fig.add_subplot(111)
		M = self.IG.getMicroenv(i)

		xgrid = M[0, :].reshape(self.IG.numy,self.IG.numx)
		ygrid = M[1, :].reshape(self.IG.numy,self.IG.numx)

		N = self.IG.getSubstrate(i, substrate_name)
		if N.shape[0] != 1:
			raise ValueError("Plotter.animate: Expected 2D microenvironment. Check input files.")

		# Actualy plot the surface
		fig.clear()
		ax_sim = fig.add_subplot(111)

		# Create Titles for plots
		ax_sim.set_title("PhysiCell: Substrate " + substrate_name + ":  $\\rho(x,y)$", fontsize=15)

		# Fill subplots with information
		# ax_sim.plot_surface(xgrid, ygrid, N[0], cmap=self.cmap)
		x = xgrid[0]
		y = ygrid[:,0]
		N_ticks_x = min(len(x), 6 + len(x) % 2)
		N_ticks_y = min(len(y), 6 + len(x) % 2)
		ticks_step_x = int(len(x)/N_ticks_x)
		ticks_step_y = int(len(y)/N_ticks_y)
		ax_sim.set_xticks(ticks_step_x*np.arange(N_ticks_x+1), labels=[x[min(len(x), ticks_step_x*i)] for i in range(N_ticks_x+1)])
		ax_sim.set_yticks(ticks_step_y*np.arange(N_ticks_y+1), labels=[y[min(len(x), ticks_step_y*i)] for i in range(N_ticks_y+1)][::-1])
		ax_sim.set_xlabel("x [micron]")
		ax_sim.set_ylabel("y [micron]")
		ax_sim.xaxis.set_label_position('top')
		ax_sim.xaxis.tick_top()
		pos = ax_sim.imshow(np.flip(N[0], 0), interpolation='bicubic', cmap=self.cmap, vmin=self.N_min, vmax=self.N_max)
		cbar = fig.colorbar(pos, ax=ax_sim, format='%E')#, cax=fig.add_axes([0.1, (mi-self.N_min)/(self.N_max-self.N_min), 0.2, (ma-mi)/(self.N_max-self.N_min)]))
		cbar.minorticks_on()

		fig.tight_layout()
		fig.savefig(folder + "/" + "/plot_substr_" + substrate_name + "_" + "%08d" % i + ".png")
		plt.close(fig)


	def exportPictures(self, substrate_name, folder="./pictures"):
		'''Creates Pictures for each simulation step and saves them to specified folder.'''
		# Define maximum and minimum value of plotting range for whole plotting time
		N_all = self.IG.getSubstrateSeries(substrate_name)
		self.N_max = max(np.array(N_all).flatten())
		self.N_min = min(np.array(N_all).flatten())

		# Initialize pool of workers
		pool = mp.Pool(mp.cpu_count()-8)

		partial_func = partial(self.exportSinglePicture,substrate_name, folder)

		with pool as p:
			r = list(tqdm.tqdm(p.imap(partial_func , range(len(N_all))), total=len(N_all)))