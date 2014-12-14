import sys
import subprocess
#sys.popen("./project2 data/grid3.nodes data/grid3.tets")i
avg = 0.0

arr = []
for a in sys.argv[1:-1]:
  arr.append(a)

N = int(sys.argv[-1])

for i in range(N):
	p = subprocess.Popen(arr,stdout=subprocess.PIPE)
	out,err = p.communicate()
	avg += float(out[2:-2])

print avg/N
