import numpy
import scipy
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Straighten raw 3D image data of worms.')
parser.add_argument('filename', type=str);
parser.add_argument('--width',  type=int, required=True);
parser.add_argument('--height', type=int, required=True);
args = parser.parse_args()
width = args.width
height = args.height
num_slices = os.path.getsize(args.filename)/2/width/height

data = numpy.memmap(args.filename, dtype=numpy.uint16, mode='r')
#step 0 done

sampled_data = data[numpy.random.randint(data.size, size=5000)]
std = numpy.std(sampled_data)
mean = numpy.mean(sampled_data)
threshold = mean+std
upper_threshold = threshold+std
print "The threshold is: ", threshold
#step 1 done

data.shape=(num_slices, args.height, args.width)
w = numpy.empty((0, 4), dtype=numpy.uint32)
while(w.shape[0] < 150):
  ind = (numpy.random.randint(data.shape[0]),numpy.random.randint(data.shape[1]),numpy.random.randint(data.shape[2]))
  if(data[ind] > threshold): # and data[ind] < upper_threshold):
    w=numpy.append(w,[numpy.append(ind, [data[ind]])], axis=0)
#step 2 done

plt.imshow(data[80,:,:])
print w
plt.scatter(w[:,2], w[:,1], 12, c=w[:,0], cmap='jet')
#plt.show()

e = numpy.empty((w.shape[0], w.shape[0]))
for i in range(w.shape[0]):
  for j in range(w.shape[0]):
    e[i,j]=scipy.sqrt((w[i,0]-w[j,0])**2 + (w[i,1]-w[j,1])**2 + (w[i,2]-w[j,2])**2)
#step 3 done

#closest_node = numpy.amin(e,axis=0)

mst_v = set([numpy.random.randint(w.shape[0])])
mst_e = set()
while(len(mst_v) < w.shape[0]):
  m = numpy.inf
  m_e = (-1,-1)
  for v in mst_v:
    for v2 in range(w.shape[0]):
      if v2 not in mst_v:
        if e[v,v2] < m:
          m_e = (v,v2)
          m = e[v,v2]
  mst_e.add(m_e)
  mst_v.add(m_e[1])
#step 4 done

for e in mst_e:
  plt.arrow(w[e[0],2],w[e[0],1],w[e[1],2]-w[e[0],2],w[e[1],1]-w[e[0],1])
#plt.show()

def bfs(k,g,wt,a,o):
  if o[k] >= 0:
    o[k]=a
    for e in g:
      if e[0] == k:
        bfs(e[1],g,wt,a+wt[e[0],e[1]],o)
      if e[1] == k:
        bfs(e[0],g,wt,a+wt[e[0],e[1]],o)

bfs1=numpy.empty((w.shape[0])).fill(-1)
bfs(0,mst_e,e,0,bfs1)
bfs2=numpy.empty((w.shape[0])).fill(-1)
bfs(numpy.argmax(bfs1),mst_e,e,0,bfs2)

#def shortest_path(j,k,g,wt,o):
