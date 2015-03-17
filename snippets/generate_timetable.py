
import sys
if len(sys.argv) != 4:
	print "usage: loadsteps = sys.argv[1], max_disp = sys.argv[2], filename = sys.argv[3]\n"
	exit()
	
loadsteps= int(sys.argv[1])
max_disp = float(sys.argv[2])
filename = sys.argv[3]

f = open(r"%s.time"%filename, "w")
f.write('step \t time\n')
for i in range(loadsteps+1):
	time = float(i)/loadsteps*max_disp
	f.write('%d \t %f\n'%(i,time))
f.close()
