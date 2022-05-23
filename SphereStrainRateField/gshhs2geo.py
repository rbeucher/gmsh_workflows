import struct
import array
import math
import sys

######################################################
# # parameters
# #####################################################

## input file

# http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.2.0/gshhs+wdbii_2.2.0.tbz
# use gshhs_[cilhf].b for different resolution input file
fgshhs = open("gshhs/gshhs_i.b", "rb")

## define the resolution 

# input lat lon in radian, output in meter
def edgeLength(lat, lon):
  ## uniform 100 km
  #return 1e5
  ## 10km in the meditarrean sea
  lat = lat * 180 / math.pi
  lon = lon * 180 / math.pi
  #if (lon > 350 or lon < 40) and (lat > 20 and lat < 50) :
  return 1e4
  #else :
  #  return -1

## output parameters

# geo file
fgeo = open("world.geo", "w")

# svg file
fsvg = open("world.svg", "w")

# svg output coordinates (input in cartesian coordinates on a sphere of radius 1)
# output range should be [0, 1000] x [0, 1000]
def svgOutputCoordinates(x, y, z):
  ## stereographic output
  #u = -x / (1 + z)
  #v = -y / (1 + z)
  #return (u + 10) * 50, (v + 10) * 50 
  ## longitude-latitude output
  lat = math.asin(z)
  lon = math.atan2(y, x)
  f = 1000 / (2 * math.pi)
  return (lon + math.pi) * f, (math.pi / 2 - lat) * f

# line width in the svg file
strokeWidth = 0.05

######################################################
# # end parameters
# #####################################################

try :
  sys.path.append("ann")
  import ann
  import numpy as np
  haveAnn = True
except:
  print("WARNING : ann wrapper not found, it will be slow for > 10k nodes")
  haveAnn = False

R = 6.371e6;

def xyz2latlon(x, y, z):
  lat = math.asin(z)
  lon = math.atan2(y, x)
  return lat, lon

def latlon2xyz(lat, lon):
  x = math.cos(lat) * math.cos(lon)
  y = math.cos(lat) * math.sin(lon)
  z = math.sin(lat)
  return x, y, z

binh = fgshhs.read(11 * 4)
allnodes=[]
allloops=[]

# read GSHHS data
while binh:
  id, n, source,  cross_greenwich, version, t,  west, east, south, north, area, area_full, container, ancestor = struct.unpack('>iiBBBBiiiiiiii', binh)
  if (version != 9) :
    print('Wrong GSHHS file version (2.2.0 required)')
    exit()
  vv = array.array('i')
  vv.fromfile(fgshhs, n * 2)
  if (t == 1) :
    vv.byteswap()
    toRad = 1e-6 * math.pi / 180
    reset = True
    split = False
    for i in range(n):
      if (reset):
        lastxyz = [0, 0, 0]
        lastlc = 0
        nodesXYZ = []
        reset = False
      lat = vv[2 * i + 1] * toRad 
      lon = vv[2 * i + 0] * toRad
      xyz = latlon2xyz(lat, lon)
      lc = edgeLength(lat, lon);
      dd = math.sqrt((xyz[0]-lastxyz[0])**2 + (xyz[1]-lastxyz[1])**2 + (xyz[2]-lastxyz[2])**2)
      if lc > 0 and (i == 0 or  dd > (lc + lastlc) / (2 * R)):
        nodesXYZ.append(xyz)
        lastxyz = xyz
        lastlc = lc
      end = False
      if (lc < 0) :
        split = True
      if (lc < 0 or i == n - 1) :
        if len(nodesXYZ) > 5 :
          allloops.append([len(allnodes), len(allnodes) + len(nodesXYZ), not split])
          allnodes.extend(nodesXYZ)
        if len(nodesXYZ) > 0:
          reset = True
  binh = fgshhs.read(11 * 4)

deleted = [0] * len(allnodes)
    
# delete close points
if haveAnn :
  k = 10
  tree = ann.kdtree(np.array(allnodes))
  points = np.zeros((k), dtype = np.int32)
  dist = np.zeros((k))
  for i in range(len(allnodes)):
    p = allnodes[i]
    tree.search(np.array(p), points, dist)
    lat, lon = xyz2latlon(p[0], p[1], p[2])
    lci = edgeLength(lat, lon)
    for j in range(1,k):
      x, y, z = allnodes[points[j]]
      lat, lon = xyz2latlon(x, y, z)
      dd = (0.8 * (edgeLength(lat, lon) + lci) / (2 * R)) ** 2
      if dist[j] < dd :
        deleted[i] = 1
        deleted[points[j]] = 1
else :
  for i in range (len(allnodes)) :
    xi = allnodes[i]
    lat, lon = xyz2latlon(xi[0], xi[1], xi[2])
    lci = edgeLength(lat, lon)
    for j in range(i):
      xj = allnodes[j]
      lat, lon = xyz2latlon(xj[0], xj[1], xj[2])
      dd = (0.8 * (edgeLength(lat, lon) + lci) / (2 * R)) ** 2
      dx = [xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]]
      if (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] < dd) :
        deleted[i] = 1
        deleted[j] = 1

# output geo and svg file

fgeo.write("""IP = newp;
IL = newl;
IS = news;
IF = newf;
Point(IP + 0) = {0, 0, 0};
Point(IP + 1) = {0, 0, %e};
PolarSphere(IS) = {IP + 0, IP + 1};
""" % R)
fsvg.write("""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="1000" height="1000">
""")
ip = 2
il = 0
rp = 0
for loop in allloops :
  realLength = 0
  for p in range(loop[0], loop[1]) :
    if deleted[p] == 0 :
      realLength = realLength + 1
  if realLength < 5 :
    continue
  firstp = ip
  fsvg.write("<path\n d=\"M ")
  for p in range(loop[0], loop[1]):
    rp = rp + 1
    if deleted[p] != 0 :
      continue
    [x, y, z] = allnodes[p]
    #stereographic coordinates
    u = -x / (1 + z)
    v = -y / (1 + z)
    fgeo.write("Point(IP + %i) = {%.16e, %.16e, %.16e}; // %i\n" % (ip, u, v, 0., rp -1))
    u, v = svgOutputCoordinates(x, y, z)
    fsvg.write(" %.5f,%.5f " % (u, v))
    ip = ip + 1
  if (loop[2]):
    fgeo.write("BSpline(IL + %i) = {IP + %i : IP + %i, IP + %i};\n" % (il, firstp, ip - 1, firstp))
    fsvg.write(" z\"\n")
  else :
    fgeo.write("BSpline(IL + %i) = {IP + %i : IP + %i};\n" % (il, firstp, ip - 1))
    fsvg.write("\"\n")
  fsvg.write("style=\"fill:none;stroke:#000000;stroke-width:%f;\"/>\n" % strokeWidth)
  il = il + 1
fgeo.write("Field[IF] = Attractor;\nField[IF].NodesList = {IP + 2 : IP + %i};\n" % (ip - 1));
fsvg.write("</svg>\n");
