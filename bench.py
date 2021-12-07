#!/bin/python3
import xml.etree.ElementTree as ET
import os

# Name of the config file to write the following parameters into
xml_file = "config/PhysiCell_settings.xml"
output_dir = "./logs_V"
project_name = "secretion_project"

# Define parameters to benchmark against
thread_configs = [1,2,4,8,12,16]
cell_configs = [25, 36, 49, 64, 81, 100, 121, 144, 225, 400, 625, 900, 1225, 1600, 2025, 2500]

# Determines the number of runs over which to average
N_runs = 2
total_runs = N_runs*len(thread_configs)*len(cell_configs)

# For printing with indents
pr_buff = 12
tabs = 3
# Counts the number of processes finished
i = 1

# The final line of the project gives us the time it took to complete.
# This can be parsed by this simple function.
def parseOutputLine(string):
	# 0 days, 0 hours, 0 minutes, and 0.0000 seconds
	es = string.split(" ")
	result = float(es[7]) + 60*(float(es[4]) + 60*(float(es[2]) + 24*float(es[0])))
	return result

def setParameters(xml_f, n_cells, num_threads):
	# Get the whole tree and root of the xml file
	tree = ET.parse(xml_f)
	root = tree.getroot()

	# Search in xml file and set number of cells
	parameters = root.find("user_parameters")
	N_cells = parameters.find("number_of_cells")
	N_cells.text = str(n_cells)
	# Again for threads
	parallel = root.find("parallel")
	omp_num_threads = parallel.find("omp_num_threads")
	omp_num_threads.text = str(num_threads)

	# Write to the file
	tree.write(xml_f)

def benchmarkProject(runs, filename):
	for j in range(runs):
		print(2*tabs*" " + "[" + str(i) + "/" + str(total_runs) + "]" + (pr_buff-2*tabs-len(str(total_runs))-len(str(i)))*" " + " Writing to file " + filename)
		os.system("./" + project_name + " >> " + filename)

		# Open the corresponding logfile (only reading)
		file = open(filename,"r")
		# Times that the project took to complete
		times = []
		# Works as follows: Mark that last line was found when regex below is matching
		# In the next iteration (line) the time will be given and parsed
		foundit = False
		for line in file:
			if foundit == True:
				times.append(parseOutputLine(line))
				foundit = False
			if "Total simulation runtime: " in line:
				foundit = True
		i += 1
	return times

def writeSummary(times, threads, cells):
	summary.write("Results for Configuration:\n")
	summary.write("Threads: " + str(threads) + "\n")
	summary.write("Cells:   " + str(cells) + "\n")

	# Average over the different runs (for the same configuration)
	# and store results in summary and matrix file
	times_average = sum(times)/len(times)
	for l, time in enumerate(times):
		summary.write("Run " + str(l) + ": " + str(time) +"\n")
	summary.write("Average: " + str(times_average) + "\n\n")
	matrix.write(str(times_average) + " ")

def checkIfDirComplete(dir):
	# If the dir exists check if every file for each cell-thread config is present
	if os.path.isdir(dir):
		print("Testing dir " + dir + " for completeness")
		comp = True
		for cells in cell_configs:
			for threads in thread_configs:
				filename = dir + "logs-T" + str(threads) + "-C" + str(cells) + ".txt"
				# This will yield a "False" if one file is missing.
				# Then return False afterwards
				comp *= os.path.isfile(filename)
		print("The directory "+ dir + " was " + comp*"complete" + (1-comp)*"incomplete")
		return comp
	# If it does not exist, create it and return false to start calculations
	else:
		print("Creating directory " + dir)
		# os.mkdir(dir)
		return False

c_j = 0
while checkIfDirComplete(output_dir + str(c_j) + "/"):
	c_j += 1
output_dir = output_dir + str(c_j) + "/"

# Stores the averages results in a matrix file and total logs in a logfile
matrix = open(output_dir + "log_matrix.txt","a")
summary = open(output_dir + "logs.txt","a")

# Iterate over all possible parameter combinations
for m, cells in enumerate(cell_configs):
	print("[" + str(m+1) + "/" + str(len(cell_configs)) + "]" + (pr_buff-len(str(len(cell_configs)))-len(str(m)))*" " + " Using " + str(cells) + "/" + str(cell_configs[-1]) + " cells")
	for n, threads in enumerate(thread_configs):
		filename = output_dir + "logs-T" + str(threads) + "-C" + str(cells) + ".txt"
		# Check if file exists and if so do nothing
		# This assumes previous calculations already created the file and need not be computed again
		if os.path.isfile(filename):
			existed = True
			print("Already finished job " + str(i) + "/" + str(N_runs*len(thread_configs)*len(cell_configs)))
			i += N_runs
		else:
			existed = False
			# Sets the parameters for the next benchmark runs
			setParameters(xml_file, cells, threads)

			# Display message for current run
			print(tabs*" " + "[" + str(n+1) + "/" + str(len(thread_configs)) + "]" + (pr_buff-tabs-len(str(len(thread_configs)))-len(str(n)))*" " + " Using " + str(threads) + "/" + str(thread_configs[-1]) + " threads")

			# Also display messages for every run over which is averaged
			times = benchmarkProject(N_runs, filename)

			# Write to summary file
			writeSummary(times, threads, cells)
	# Only add a new line if the logfile of the last configuration did previously not exist.
	# Otherwise we will continue the line
	if existed == False:
		matrix.write(str("\n"))
