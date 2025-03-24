'''
This notebook is to analytically show the overhead of message aggregation
protocols used in our \kmer counting workload.
'''
from defaultplot import *
from params import *
from math import *
from matplotlib import pyplot as plt
import numpy as np

# Output of each function is memory used in Bytes
def kmerbytes(k):
  return (2**(ceil(log2(2 * k)))) / 8

def l01d(p, c0):
  return (2 * c0 * p)

def l02d(p, c0):
  return (2 * c0 * p**(0.5))

def l03d(p, c0):
  return (2 * c0 * p**(1/3))

def l1(c1, pkt_load):
  return (c1 * pkt_load)

def l2(c2, p, k):
  pkt = (c2 * kmerbytes(k)) + 8 # last 8 bytes for the metadata per packet
  return (p * pkt)

def l3(c3, k):
  return (c3 * kmerbytes(k))

def total1d(p, c0, c1, c2, c3, k):
  mem3 = l3(c3, k)
  mem2 = l2(c2, p, k)
  pktsize = mem2 / p
  mem1 = l1(c1, pktsize)
  mem0 = l01d(p, c0)
  return (2 * mem0) + (2 * mem2) + mem1 + mem3

def total2d(p, c0, c1, c2, c3, k):
  mem3 = l3(c3, k)
  mem2 = l2(c2, p, k)
  pktsize = mem2 / p
  mem1 = l1(c1, pktsize)
  mem0 = l02d(p, c0)
  return (2 * mem0) + (2 * mem2) + mem1 + mem3

def total3d(p, c0, c1, c2, c3, k):
  mem3 = l3(c3, k)
  mem2 = l2(c2, p, k)
  pktsize = mem2 / p
  mem1 = l1(c1, pktsize)
  mem0 = l03d(p, c0)
  return (2 * mem0) + (2 * mem2) + mem1 + mem3

def kmermem(reads, readlen, k, p, cov):
  input_bytes = reads * readlen
  output_bytes = reads * (readlen - k + 1) * kmerbytes(k)
  return (input_bytes + output_bytes) / p

x = []
ykmr = []

yl1mem = []
yl2mem = []
yl3mem = []

yl01dmem = []
yl02dmem = []
yl03dmem = []

y1dtotal = []
y2dtotal = []
y3dtotal = []

while (p <= 2**14):
  # All memory is in bytes

  # Memory required for the algorithm itself
  kmr = kmermem(reads, readlen, k, p, cov) / (1024 * 1024)

  # Memory required for l1 -> l3 layers of aggregation
  l3mem = l3(c3, k)     / (1024 * 1024)
  l2mem = 2 * l2(c2, p, k)  / (1024 * 1024)

  pktsize = l2mem / (2 * p) # already in KB
  l1mem = l1(c1, pktsize)

  # Memory required for different conveyor aggregation
  l01dmem = 2 * l01d(p, c0) / (1024 * 1024)
  l02dmem = 2 * l02d(p, c0) / (1024 * 1024)
  l03dmem = 2 * l03d(p, c0) / (1024 * 1024)

  # Get the total Memory required
  y1dmem = total1d(p, c0, c1, c2, c3, k) / (1024 * 1024)
  y2dmem = total2d(p, c0, c1, c2, c3, k) / (1024 * 1024)
  y3dmem = total3d(p, c0, c1, c2, c3, k) / (1024 * 1024)

  # Time to perform. some assert statements
  assert np.allclose(y1dmem, l01dmem + l1mem + l2mem + l3mem)
  assert np.allclose(y2dmem, l02dmem + l1mem + l2mem + l3mem)
  assert np.allclose(y3dmem, l03dmem + l1mem + l2mem + l3mem)

  # X-axis is the number of processors
  x.append(p)

  # The different Y-axis information
  yl1mem.append(l1mem)
  yl2mem.append(l2mem)
  yl3mem.append(l3mem)
  yl01dmem.append(l01dmem)
  yl02dmem.append(l02dmem)
  yl03dmem.append(l03dmem)
  ykmr.append(kmr)

  y1dtotal.append(y1dmem)
  y2dtotal.append(y2dmem)
  y3dtotal.append(y3dmem)

  # Double the processor count
  p *= 2

# make everything np.array for simplicity
x        = np.array(x)
ykmr     = np.array(ykmr)
yl1mem   = np.array(yl1mem)
yl2mem   = np.array(yl2mem)
yl3mem   = np.array(yl3mem)
yl01dmem = np.array(yl01dmem)
yl02dmem = np.array(yl02dmem)
yl03dmem = np.array(yl03dmem)
y1dtotal = np.array(y1dtotal)
y2dtotal = np.array(y2dtotal)
y3dtotal = np.array(y3dtotal)

# Time to plot the data now
x_pos = np.arange(len(x))

# Decide the colors
color0      = '#FF8C00'    # DarkOrange
color1      = '#3CB371'    # MediumSeaGreen
color2      = '#FFD700'    # Gold
color3      = '#DC143C'    # Crimson
color_algo  =  colors[8]# '#8A2BE2'    # BlueViolet
edge        = 'black'      # Black - duh !!

# Bar width
barsize = 0.2

# X-position of each column
col1pos = x_pos - 1.5 * barsize
col2pos = x_pos - 0.5 * barsize
col3pos = x_pos + 0.5 * barsize
col4pos = x_pos + 1.5 * barsize

# Make another plot
fig, ax = plt.subplots(figsize=(6, 3.5))

# minimum y-axis is 2^2
ax.set_ylim(bottom=0, top=6)

# ax.set_yscale('log', base=2)

# Make everything Giga Bytes 

ykmr_g = ykmr / (1024)
y1dtotal_g = y1dtotal / (1024)
y2dtotal_g = y2dtotal / (1024)
y3dtotal_g = y3dtotal / (1024)

# Column 1:
ax.bar(col1pos, ykmr_g, barsize, label='Algorithm',  color=color_algo, edgecolor=edge, hatch='//')
ax.bar(col1pos, y1dtotal_g, barsize, label='1D Buf', color=colors[1], edgecolor=edge,
       bottom=ykmr_g)

# Column 2:
ax.bar(col2pos, ykmr_g, barsize,  color=color_algo, edgecolor=edge, hatch='//')
ax.bar(col2pos, y2dtotal_g, barsize, label='2D Buf', color=colors[2], edgecolor=edge,
       bottom=ykmr_g)

# Column 3:
ax.bar(col3pos, ykmr_g, barsize,  color=color_algo, edgecolor=edge, hatch='//')
ax.bar(col3pos, y3dtotal_g, barsize, label='3D Buf', color=colors[0], edgecolor=edge,
       bottom=ykmr_g)

# Add labels and title
ax.set_xlabel('Processors')
ax.set_ylabel('Memory (GB)')
# ax.set_title('Total Memory Usage')
ax.set_xticks(x_pos)
ax.set_xticklabels(x)
ax.legend()

# Save the plot
plt.tight_layout()
plt.savefig('figures/total_mem.svg', format='svg')