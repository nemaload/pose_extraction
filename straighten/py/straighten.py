import numpy
import numpy.linalg
import ctypes
import scipy
import scipy.interpolate
import scipy.ndimage
import scipy.weave
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
while(w.shape[0] < 80):
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
beta = 0
gamma = 7.0

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
  return (sss/new.shape[0])

ps = numpy.empty((0, 4), dtype=numpy.uint32)
while(ps.shape[0] < 800):
  ind = (numpy.random.randint(data.shape[0]),numpy.random.randint(data.shape[1]),numpy.random.randint(data.shape[2]))
  if(data[ind] > mean): # and data[ind] < upper_threshold):
    ps=numpy.append(ps,[numpy.append(ind, [data[ind]])], axis=0)

iteration_count = 0
while(ens_d(ck_new,ck) > 40):
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
    if(k>1):
      nk[k-1,0]+=p[3]
      ipp[k-1]+=p[3]*p[0:3]
    nk[k,0]+=p[3]
    ipp[k]+=p[3]*p[0:3]
    if(k<nk.shape[0]-1):
      nk[k+1,0]+=p[3]
      ipp[k+1]+=p[3]*p[0:3]

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

#Step 9 - mmaping output image

new_data = numpy.memmap(args.filename + ".small", dtype=numpy.uint16, mode='write',shape=(data.shape[0]/20,data.shape[1]/20,data.shape[2]/10))

#Step 8 - modifying the image based on this backbone

((tck,u),fp,err,msg) = scipy.interpolate.splprep(ck_new.transpose(),full_output=1)
print "Beginning step 8/9..."

def splinterp(output_coords):
  if(output_coords[2] == 0):
    print "Computing point", output_coords
  x = output_coords[2]
  bp = scipy.interpolate.splev(x,tck)
  bd = scipy.interpolate.splev(x,tck,der=1)
  znv = [1,0,0]
  ynv = [0,1,0]
  zdv = numpy.cross(bd,ynv)
  ydv = numpy.cross(bd,znv)
  zndv = zdv/numpy.linalg.norm(zdv)
  yndv = ydv/numpy.linalg.norm(ydv)
  if(zndv[0] < 0):
    zndv*=-1
  return bp + zndv*(output_coords[0]-data.shape[0]/2) + yndv*(output_coords[1]-data.shape[1]/2)

#scipy.ndimage.geometric_transform(data,splinterp,output=new_data,order=0,prefilter=False)

n1k, n1j, n1i = data.shape
n2k, n2j, n2i = new_data.shape
cdata = data.ctypes.data
cnew_data = new_data.ctypes.data
cdata_high = cdata / 2**32
cnew_data_high = cnew_data / 2**32
cdata=cdata % 2**32
cnew_data = cnew_data % 2**32
for oi in range(n2i):
  print "Computing spline point", oi
  bp = scipy.interpolate.splev(oi*1.0/new_data.shape[2],tck)
  bd = scipy.interpolate.splev(oi*1.0/new_data.shape[2],tck,der=1)
  znv = [1,0,0]
  ynv = [0,1,0]
  zdv = numpy.cross(bd,ynv)
  ydv = numpy.cross(bd,znv)
  zndv = zdv/numpy.linalg.norm(zdv)
  yndv = ydv/numpy.linalg.norm(ydv)
  if(zndv[0] < 0):
    zndv*=-1
  print cnew_data
  code =  """
          #line 330 "straighten2.py"
          unsigned long long pnew_data = (unsigned long long)cnew_data;
          unsigned short* new_data = (unsigned short*)pnew_data;
          unsigned short* data = (unsigned short*)cdata;
          for(int ok=0; ok<n2k; ok++) {
            printf("Row %d\\n", ok);
            for(int oj=0; oj<n2j; oj++) {
              int iip[3];
              int in_bounds = 1;
              for(int ov=0; ov<3; ov++) {
                printf("ov=%d\\n",ov);
                iip[ov]=(int)((double)bp[ov]+zndv[ov]*(ok*(n2k*1.0/n1k)-n1k/2.0)+yndv[ov]*(oj*(n2j*1.0/n1j)-n1j/2.0));
              }
              if(iip[0] >= 0 && iip[1] >= 0 && iip[2] >= 0 && iip[0] < n1k && iip[1] < n1j && iip[2] < n1i) {
                printf("Trying (%d,%d,%d)\\n",iip[0],iip[1],iip[2]);
                //new_data(ok,oj,oi) = data(iip[0],iip[1],iip[2]);
              } else {
                printf("Zeroing (%d,%d,%d)\\n",ok,oj,oi);
                printf("%llu[%d]\\n",new_data,ok*n2j*n2i+oj*n2i+oi);
                new_data[ok*n2j*n2i+oj*n2i+oi] = 0;
              }
            }
          }
          """
  scipy.weave.inline(code,['oi','bp','zndv','yndv','n1i','n1j','n1k','n2i','n2j','n2k','cnew_data','cnew_data_high','cdata','cdata_high'])
  #for oj in range(new_data.shape[1]):
    #print "Row", oj
    #for ok in range(new_data.shape[0]):
      #ip = bp + zndv*(ok*20-data.shape[0]/2) + yndv*(oj*20-data.shape[1]/2)
      #(ik,ij,ii) = numpy.floor(ip).astype(int)
      #if(ii >= 0 and ij >= 0 and ik >= 0 and ik < data.shape[0] and ij < data.shape[1] and ii < data.shape[2]):
        #new_data[ok,oj,oi]=data[ik,ij,ii]
      #else:
        #new_data[ok,oj,oi]=0

plt.imshow(new_data[4,:,:])

plt.show()
