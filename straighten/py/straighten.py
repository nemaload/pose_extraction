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
print "Step 0 done!"

sampled_data = data[numpy.random.randint(data.size, size=5000)]
std = numpy.std(sampled_data)
mean = numpy.mean(sampled_data)
threshold = mean+std
upper_threshold = threshold+std
#print "The threshold is: ", threshold
print "Step 1 done!"

data.shape=(num_slices, args.height, args.width)
w = numpy.empty((0, 4), dtype=numpy.uint32)
while(w.shape[0] < 250):
  ind = (numpy.random.randint(data.shape[0]),numpy.random.randint(data.shape[1]),numpy.random.randint(data.shape[2]))
  if(data[ind] > threshold): # and data[ind] < upper_threshold):
    w=numpy.append(w,[numpy.append(ind, [data[ind]])], axis=0)
print "Step 2 done!"

plt.imshow(data[80,:,:])
#print w
plt.scatter(w[:,2], w[:,1], 12, c=w[:,0], cmap='jet')
#plt.show()

ew = numpy.empty((w.shape[0], w.shape[0]))
for i in range(w.shape[0]):
  for j in range(w.shape[0]):
    ew[i,j]=scipy.sqrt((w[i,0]-w[j,0])**2 + (w[i,1]-w[j,1])**2 + (w[i,2]-w[j,2])**2)
print "Step 3 done!"

#closest_node = numpy.amin(e,axis=0)

mst_v = set([numpy.random.randint(w.shape[0])])
mst_e = set()
while(len(mst_v) < w.shape[0]):
  m = numpy.inf
  m_e = (-1,-1)
  for v in mst_v:
    for v2 in range(w.shape[0]):
      if v2 not in mst_v:
        if ew[v,v2] < m:
          m_e = (v,v2)
          m = ew[v,v2]
  mst_e.add(m_e)
  mst_v.add(m_e[1])
print "Step 4 done!"

def arrow(u,v):
  plt.arrow(w[u,2],w[u,1],w[v,2]-w[u,2],w[v,1]-w[u,1])

#for e in mst_e:
#  arrow(e[0],e[1])
#plt.show()

def bfs(k,g,wt,a,o):
  if o[k] < 0:
    o[k]=a
    for e in g:
      if e[0] == k:
        bfs(e[1],g,wt,a+wt[e[0],e[1]],o)
      if e[1] == k:
        bfs(e[0],g,wt,a+wt[e[0],e[1]],o)

bfs1=numpy.empty((w.shape[0]))
bfs1.fill(-1.0)
bfs(0,mst_e,ew,0,bfs1)
bfs2=numpy.empty((w.shape[0]))
bfs2.fill(-1)
bfs(numpy.argmax(bfs1),mst_e,ew,0,bfs2)

tip1=numpy.argmax(bfs1)
tip2=numpy.argmax(bfs2)

dist = dict()
for v in range(w.shape[0]):
  dist[v]=numpy.Inf

previous = dict()
dist[tip1]=0
q = mst_v.copy()
while len(q) > 0:
  u = -1
  dist_u = numpy.Inf
  for v in q:
    if dist[v] < dist_u:
      u = v
      dist_u = dist[v]
  if(u==-1):
    break
  q.remove(u)
  for (x,y) in mst_e:
    if x == u:
      if y in q:
        alt = dist[x] + ew[x,y]
        if alt < dist[y]:
          dist[y] = alt
          previous[y] = x
    if y == u:
      if x in q:
        alt = dist[y] + ew[y,x]
        if alt < dist[x]:
          dist[x] = alt
          previous[x] = y

print "Step 5 done!"

it = tip2
while(it in previous):
  arrow(it,previous[it])
  it = previous[it]

#plt.show()
