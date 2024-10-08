from experiments import *
from params import *
from kcount import *
from defaultplot import *

# Perform strong scaling experiments on the S32 dataset and observe
# trend of hardware resource utilization 
num_reads32 = 1431655750
P = 4
start = int(math.log2(P))

intra = []
inter = []
comp = []

total = []
procs = []

i = 0
while P < 2**16:
  cmpu = cpu_compute(m, num_reads32, k, P, L_cpu, Z_cpu, cnode_cpu)
  itra = cpu_intranode(m, num_reads32, k, P, L_cpu, Z_cpu, bmem_cpu, A)
  interr = cpu_internode(m, num_reads32, k, P, blink)
  
  assert(cmpu > 0)
  assert(itra > 0)
  assert(interr > 0)

  tot = cmpu + itra + interr
  
  comp.append(cmpu/tot)
  intra.append(itra/tot)
  inter.append(interr/tot)
  total.append(tot)

  procs.append(P)

  P *= 2

# Decide the colors
color0      = '#FF8C00'    # DarkOrange
color1      = '#3CB371'    # MediumSeaGreen
color2      = '#FFD700'    # Gold
color3      = '#DC143C'    # Crimson
color4      = '#8A2BE2'    # BlueViolet
edge        = 'black'      # Black - duh !!

# Plotting the marks + line plot
ind = np.arange(len(procs))  # the x locations for the groups

fig, ax = plt.subplots(figsize=(5, 3))

# make the y-axis log 2 scale 

# Create a stacked bar plot
width = 0.35  # the width of the bars

p1 = ax.bar(ind, intra, width, label='Intra-node', color=color0, edgecolor=edge)
p2 = ax.bar(ind, inter, width, bottom=intra, label='Inter-node', color=color4, edgecolor=edge)
p3 = ax.bar(ind, comp, width, bottom=[i+j for i,j in zip(intra, inter)], label='Compute', color=color1, edgecolor=edge)

# make sure y-axis has 0.5 on it
ax.set_ylim([0, 1])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])

# make a horizontal black line y = 0.5, and mark it 
ax.axhline(y=0.5, color='black', linewidth=0.5)

ax.set_xlabel('Nodes')
ax.set_ylabel('Predicted Time (s)')
ax.set_title('Synthetic 32', fontsize='large')
# ax.set_yscale('log', base=2)
# yticks = [2**i for i in range(7, 10)]
# ax.set_yticks(yticks)
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks(ind)
ax.set_xticklabels([f'$2^{{{i+start}}}$' for i in range(len(procs))], fontsize='small')
ax.tick_params(axis='y', labelsize='small')
ax.legend()

plt.savefig('figures/sscale32model.svg', format='svg')
# plt.show()