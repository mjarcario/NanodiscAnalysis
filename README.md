# NanodiscAnalysis
A set of procedures to analyze molecular dynamics trajectories of protein-embedded nanodisc assemblies.

This repository contains a set of Tcl and Python scripts to analyze the properties of incorporated proteins 
and lipids in a nanodisc assembly. These scripts have been employed in VMD (1.9.4a55 using Tcl version 8.5.6) and 
Python 3 (versions 3.6.8 and 3.11.4). The main analysis script, which contains the majority of procedures used to analyze
the nanodisc assemblies is contained in "nanodiscAnalysis.tcl"

"procedures.tcl" contains generic procedures for the analysis of molecular dynamics trajectories. This scripts is called 
in the main body of "nanodiscAnalysis.tcl". Ensure that the relative path for "procedures.tcl" is correct within "nanodiscAnalysis.tcl".

The "getNanodiscEllipse" function in nanodiscAnalysis.tcl requires an external python script called "ellipse.py"
which takes the coordinates of the atoms and produces an ellipse of best-fit along with the major and minor axes.
