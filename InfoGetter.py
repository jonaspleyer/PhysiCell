#!/bin/python3

import xml.etree.ElementTree as ET
import os
import scipy.io

class InfoGetter:
	def __init__(self, folder=os.getcwd() + "/output"):
		# Store folder name of output files
		self.output_folder_full = folder

		# Check if the supplied output self.output_folder_full exists
		if not os.path.isdir(self.output_folder_full):
			raise RuntimeError("InfoGetter.__init__: Folder " + self.output_folder_full + " does not appear to exist")

		# Store the current working directory
		self.cwd = os.getcwd()

		# Store Information about where files are presented
		self.getFilesInfo()

		# Store information about the coordinates and sizes
		self.getDomainInfo()

	def getFilesInfo(self):
		'''Store Information about where files are presented'''
		self.files_init = [filename for filename in os.listdir(self.output_folder_full) if "init" in filename]
		self.files_all = os.listdir(self.output_folder_full)
		self.files_all_reduced = [f for f in self.files_all if "final" not in f and "init" not in f]

		# Store filenames of xml files
		self.files_xml_all = sorted([f for f in self.files_all if ".xml" in f])
		self.files_xml_reduced = [f for f in self.files_xml_all if "final" not in f and "init" not in f]

		# Store filenames of matlab files
		self.files_mat_all = sorted([f for f in self.files_all if ".mat" in f])
		self.files_mat_micros = sorted([f for f in self.files_all_reduced if "microenvironment" in f and ".mat" in f])
		self.files_mat_cells_physicell = sorted([f for f in self.files_all_reduced if "cells_physicell" in f and ".mat" in f])
		self.files_mat_cells = sorted([f for f in self.files_all_reduced if "cells" in f and "physicell" not in f and ".mat" in f])

	def getDomainInfo(self):
		''''Stores information about domain size of the simulation. Utilizes the initial.xml file and assumes that the domain does not change.'''
		X_Node = self.getXML_Node("initial","microenvironment/domain/mesh/x_coordinates")
		Y_Node = self.getXML_Node("initial","microenvironment/domain/mesh/y_coordinates")
		Z_Node = self.getXML_Node("initial","microenvironment/domain/mesh/z_coordinates")
		self.x_coordinates = X_Node.text.split(X_Node.attrib["delimiter"])
		self.y_coordinates = Y_Node.text.split(Y_Node.attrib["delimiter"])
		self.z_coordinates = Z_Node.text.split(Z_Node.attrib["delimiter"])
		self.numx = len(self.x_coordinates)
		self.numy = len(self.y_coordinates)
		self.numz = len(self.z_coordinates)

	def getFileName(self, i, name_list):
		'''Determines if given input is acceptable and then returns the name of the file for the specified simulation step.'''
		if type(i) == str and ("final" == i or "initial" == i):
			results = [n for n in name_list if i in n]
			if len(results)>1:
				raise RuntimeError("InfoGetter.getFileName: Multiple files for Index " + str(i) + " were found. Check your input files.")
			else:
				return self.output_folder_full + "/" + results[0]
		elif type(i) == int and i>-len(name_list) and i<len(name_list):
			return self.output_folder_full + "/" + name_list[i]
		elif type(i) == int:
			raise IndexError("InfoGetter.getFileName: Index " + str(i) + " out of bounds: no output file found")
		else:
			raise TypeError("InfoGetter.getFileName: Index should be integer in bounds or \"final\" or \"initial\"")

	def getXML_Node(self, i, valname):
		'''Searches the xml file of the specified simulation step for a node and returns it.'''
		xml = ET.parse(self.getFileName(i, self.files_xml_all))
		root = xml.getroot()
		node = root.find(valname)
		return node

	def getMicroenv(self, i):
		'''Get the microenvironment of the specified simulation step.'''
		name = self.getFileName(i, self.files_mat_micros)
		info_dict = {}
		scipy.io.loadmat(name, info_dict)
		M = info_dict["multiscale_microenvironment"]
		return M
	
	def getMicroenvSeries(self):
		'''Get the microenvironment series of the total simulation.'''
		return [self.getMicroenv(i) for i in range(len(self.files_mat_micros))]

	def getSubstrate(self, i, substrate_name):
		'''Returns the substrate densities in the specified simulation step.'''
		variables = self.getXML_Node("initial","microenvironment/domain/variables")
		ID = -1
		for var in variables:
			ID += 1
			if var.attrib["name"] == substrate_name:
				# ID = int(var.attrib["ID"])
				break
		if ID == -1:
			raise IndexError("InfoGetter.getSubstrateSeries: Cannot find Index for substrate " + substrate_name + " check your name in xml file")
		return self.getMicroenv(i)[ID+4,:].reshape(self.numz, self.numy, self.numx)

	def getSubstrateSeries(self, substrate_name):
		'''Generate a list of Microenvironment arrays of shape (numz, numy, numx).'''
		Subs = []
		for n in range(len(self.files_mat_micros)):
			N = self.getSubstrate(n, substrate_name)
			Subs.append(N)
		return Subs
