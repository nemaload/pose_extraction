import numpy
import scipy
import argparse
import os
import matplotlib.pyplot as plt
from pylab import *

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
while(w.shape[0] < 50):
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

def prim_arrow(p1,p2,color='black'):
  plt.arrow(p1[2],p1[1],p2[2]-p1[2],p2[1]-p1[1],color=color)

def arrow(u,v):
  prim_arrow(w[u],w[v])

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

ck = numpy.empty((0, 3))
it = tip2
while(it in previous):
  arrow(it,previous[it])
  ck=numpy.append(ck,[[w[it,0],w[it,1],w[it,2]]], axis=0)
  it = previous[it]

print "Step 6 done!"

alpha = 1
beta = -0.4
gamma = 3.0

#next 17 lines taken from Wikipedia
class Node:pass
 
def kdtree(plist, depth=0):
    if not plist:
        return
 
    axis = depth % 3
 
    plist.sort(key=lambda point: point[axis])
    median = len(plist)/2
 
    node = Node()
    node.location = plist[median]
    node.leftChild = kdtree(plist[0:median], depth+1)
    node.rightChild = kdtree(plist[median+1:], depth+1)
    return node

def distance(p1, p2):
  ss=0
  for i in range(3):
    ss+=(p1[i]-p2[i])**2
  return scipy.sqrt(ss)

def kdsearch(here, point, best=None, axis=0):
  if here == None:
    return best

  if best == None:
    best = here

  if distance(here.location,point) < distance(best.location,point):
    best = here

  left_nearer=0
  if point[axis] < here.location[axis]:
    child = here.leftChild
    left_nearer=1
  else:
    child = here.rightChild
  best = kdsearch(child,point,best,axis=(axis+1)%3)

  if abs(here.location[axis]-point[axis]) < distance(best.location,point):
    if left_nearer:
      child=here.rightChild
    else:
      child=here.leftChild
    best=kdsearch(child,point,best,axis=(axis+1)%3)

  return best

plist = []
ck_new = ck
ck = numpy.zeros_like(ck_new)

def ens_d(new, old):
  ss = 0
  sss = 0
  for i in range(new.shape[0]):
    for j in range(new.shape[1]):
      ss += (new[i,j]-old[i,j])**2
    sss += scipy.sqrt(ss)
    ss = 0
  return sss

ps = numpy.empty((0, 4), dtype=numpy.uint32)
while(ps.shape[0] < 800):
  ind = (numpy.random.randint(data.shape[0]),numpy.random.randint(data.shape[1]),numpy.random.randint(data.shape[2]))
  if(data[ind] > mean): # and data[ind] < upper_threshold):
    ps=numpy.append(ps,[numpy.append(ind, [data[ind]])], axis=0)

iteration_count = 0
while(ens_d(ck_new,ck) > 200):
  print "Iteration", iteration_count
  iteration_count+=1
  ck = ck_new

  ps = numpy.delete(ps,np.s_[0:50], 0)
  while(ps.shape[0] < 800):
    ind = (numpy.random.randint(data.shape[0]),numpy.random.randint(data.shape[1]),numpy.random.randint(data.shape[2]))
    if(data[ind] > mean): # and data[ind] < upper_threshold):
      ps=numpy.append(ps,[numpy.append(ind, [data[ind]])], axis=0)

  plist = []
  for k in range(ck.shape[0]):
    p = ck[k]
    plist.append({0: p[0], 1: p[1], 2: p[2], "k": k})
  kd = kdtree(plist)
  nk = numpy.zeros((ck.shape[0],1))
  ipp = numpy.zeros((ck.shape[0],3))

  for p in ps:
    nn = kdsearch(kd,p).location
    k = nn["k"]
    nk[k,0]+=p[3]
    ipp[k]+=p[3]*p[0:3]

  #print nk[0]
  for i in range(ck.shape[0]):
    if(nk[i,0] != 0):
      nk[i,0]=1/nk[i,0]
  #print nk[0]

  ckm1 = numpy.concatenate((ck[1:], ck[-2:-1]))
  ckm2 = numpy.concatenate((ck[2:], ck[-2:-1], ck[-2:-1]))
  ckp1 = numpy.concatenate((ck[0:1], ck[:-1]))
  ckp2 = numpy.concatenate((ck[0:1], ck[0:1], ck[:-2]))
  #print ck[-2:-1]
  #print ckm2

  #ck_new = (alpha*(1/nk)*ipp+beta*(ckm1+ckp1)+gamma*((3/2)*(ckm1+ckp1)-(1/2)*(ckm2+ckp2)))/(alpha+2*beta+2*gamma)
  ck_new = ipp
  print ck_new[0:10]
  ck_new *= nk
  print ck_new[0:10]
  ck_new *= alpha
  print ck_new[0:10]
  ck_new += beta*(ckm1+ckp1)
  print ck_new[0:10]
  ck_new += gamma*((3/2)*(ckm1+ckp1)-(1/2)*(ckm2+ckp2))
  print ck_new[0:10]
  ck_new /= (alpha+2*beta+2*gamma)
  print ck_new[0:10]

print "Step 7 done!"

plist = []
for k in range(ck_new.shape[0]):
  p = ck_new[k]
  plist.append({0: p[0], 1: p[1], 2: p[2], "k": k})
for p in ps:
  nn = kdsearch(kd,p).location
  prim_arrow(nn,p,color="green")

for k in range(ck_new.shape[0]-1):
  prim_arrow(ck_new[k],ck_new[k+1],color="red")

plt.show()
