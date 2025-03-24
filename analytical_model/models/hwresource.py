
from experiments import *
from params import *
from kcount import *
from defaultplot import *

# I now need a model that shows me the pie chart of how much time is spent 
# doing computation, intranode computation, and internode computation

# Pick the S30 Data and show what happens to different forms of hardware 
# utilization while strong scaling from 4 nodes to 128 nodes 

P = 32 # We start from 32 nodes 

comp = []
intra = []
inter = []

normalization = []

for input in inputlist:
  comp. append(cpu_compute(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu))
  intra.append(cpu_intranode(m, num_reads[input], k, P, L_cpu, Z_cpu, bmem_cpu, A))
  inter.append(cpu_internode(m, num_reads[input], k, P, blink))

for i in range(0, len(comp)):
  normalization.append(comp[i] + intra[i] + inter[i])
  comp[i] /= normalization[i]
  intra[i] /= normalization[i]
  inter[i] /= normalization[i]

# print(comp)
# print(intra)
# print(inter)

# Decide the colors
color0      = '#FF8C00'    # DarkOrange
color1      = '#3CB371'    # MediumSeaGreen
color2      = '#FFD700'    # Gold
color3      = '#DC143C'    # Crimson
color4      = '#8A2BE2'    # BlueViolet
edge        = 'black'      # Black - duh !!

# Make a pie chart for the last input
fig, ax = plt.subplots(figsize=(4, 3))

labels = ['Compute', 'Intranode', 'Internode']
sizes = [comp[-1], intra[-1], inter[-1]]
colors = [color1, color0, color4]
explode = (0, 0, 0)  # explode the first slice

ax.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
        shadow=False, startangle=90, frame=False)
ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.savefig('figures/hwpie.svg', format='svg')
# plt.show()
