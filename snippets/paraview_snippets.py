import os
import sys
from paraview.simple import *
if len(sys.argv) != 3:
	print "usage: dir = sys.argv[1], filename = sys.argv[2]\n"
	exit()
	
dir = sys.argv[1]
filename = sys.argv[2]
#dir = r"c:\Users\Dmitry\Downloads\R20H20Post"
#filename = "input"
#os.chdir(r"c:\Users\Dmitry\Downloads\R20H20Post")

vtk_ext = 'vtu'

f = open(os.path.join(dir, filename+'.time'), 'r')
f.readline()
steps = []
times = []
while (True):
	list = f.readline().split()
	if len(list) < 2:
		break
	steps.append(int(list[0]))
	times.append(float(list[1]))
f.close()

print steps
print times

r = LegacyVTKReader()
w = XMLUnstructuredGridWriter()
w.CompressorType = "ZLib"
#use base64 encoding
#w.EncodeAppendedData = 1
for i in steps:
	r.FileNames=[os.path.join(dir,'%s%d.vtk'%(filename,i))]
	w.FileName = os.path.join(dir,'%s%d.%s'%(filename, i, vtk_ext))
	w.UpdatePipeline()

#write PVD file (Paraview collection)

f = open(os.path.join(dir,'%s.pvd'%filename), 'w')
f.write('<?xml version="1.0"?>\n')
f.write('<VTKFile type="Collection" version="0.1">\n')
f.write('<Collection>\n')
for i in range(len(steps)):
	f.write('<DataSet timestep="%f" file="%s"/>\n'%(times[i], os.path.relpath(os.path.join(dir,'%s%d.%s'%(filename, steps[i], vtk_ext)) ,dir)))
f.write('</Collection>\n')
f.write('</VTKFile>\n')
f.close()
