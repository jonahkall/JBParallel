import sys
import subprocess
#sys.popen("./project2 data/grid3.nodes data/grid3.tets")i
N = 1000
avg = 0.0
for i in range(N):
	p = subprocess.Popen(["./project2", "data/grid3*"],stdout=subprocess.PIPE)
	out,err = p.communicate()
	avg += float(out[2:-2])

print avg/N
