#!/bin/python3

from argparse import ArgumentParser

# Import custom modules
from InfoGetter import *
from Plotter import *

if __name__ == '__main__':
	PL = Plotter()

	PL.exportPictures("activator")
