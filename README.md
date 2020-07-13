--------------------------
Author: Dennis J. Hazelett
Date: 7/13/2020
--------------------------

# README

Recombulator calculates the expected recombination frequency between any two loci on any chromosome of the fly genome, given as a gene name or a genomic coordinate. Recombination frequencies vary along all eukaryotic chromosomes, and these frequencies have been determined empirically with relatively high precision in the genome (Ashburner). 

## Requirements

R. Suggested: Rstudio

## Running

The easiest way to run these functions is to open "Recombulator.Rmd" in Rstudio, and run all the chunks so that the example code renders. Look at the examples and code, and change the desired input values before rerunning the chunks. Once you understand the output graphs, it should be fairly obvious what you need to do.

## To Do:

Shiny web app (help!)

# OLD README INFO

I originally wrote these functions in python as a Tkinter app, (original file date 1/27/2009). If anyone can help make these run again I would appreciate it greatly!! Here are the original README contents.

Recombulator is a suite of python functions that serve to calculate
the recombination frequency between any two loci on any chromosome of
the fly genome, given as a gene name or a genomic coordinate. The code
has been tested on a python interpreter shell on a Debian Linux with
AMD64/x86-64 in 32 bits and later in 64 bits on a dual core processor.

This directory contains two files that are necessary to run
recombulator.

recombTools.py
geneheader.txt

## REQUIREMENTS

Tkinter, python-psyco (recommended)

## RUNNING RECOMBULATOR

To run recombulator in user-friendly mode, simply execute
recombTools.py in python from the command line, and the GUI will
appear.

## NOTES 

If you have the python-psyco module installed you can uncomment the
lines pertaining to psyco at the beginning of the program.

COPYRIGHT Dannis J. Hazelett. All code described herein is distributed
under the GNU public license (GPL). A copy of the license
(gpl-3.0.txt) is included in this directory, but can also be obtained
from http://www.gnu.org/licenses/gpl.html.
