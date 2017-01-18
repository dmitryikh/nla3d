This folder contains data for solving simple 2D thermal problem with
convection and a point heat source.

The problem was taken from attached HTML files.

heat.apdl represents APDL script for Ansys Mechanical APDL to reproduce this
problem in Ansys. This script also dumps heat.cdb to be able to read mesh data
into nla3d.

temp.txt contains temperature in every node as result of steady solution.
temphistory.txt contains temperature of the upper right node (component
'PROBE' in cdb file) over the time for transient analysis.
