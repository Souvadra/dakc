from experiments import *
from params import *
from kcount import *
from defaultplot import *

# Use the standard model to predict the time taken by the algorithm in the first
# ten datasets on 8 nodes of Phoenix

P = 8 

esttime_max = []
esttime_sum = []

p1time_max = []
p1time_sum = []
p2time = []

estsortcache = []
estparsecache = []

for input in inputlist:
  esttime_max.append(cpu_time(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu, bmem_cpu, blink, A, "max"))
  esttime_sum.append(cpu_time(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu, bmem_cpu, blink, A, "sum"))

  p1time_sum.append(cpu_phase1time(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu, bmem_cpu, blink, A, "sum"))
  p1time_max.append(cpu_phase1time(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu, bmem_cpu, blink, A, "max"))
  p2time.append(cpu_phase2time(m, num_reads[input], k, P, L_cpu, Z_cpu, cnode_cpu, bmem_cpu, blink, A))

  x = N(m, num_reads[input], k, P)
  estsortcache.append(sortcomm(x, Z_cpu, L_cpu, k, A))

  estparsecache.append(t1intracachemiss(m, num_reads[input], P, L_cpu, k))

num_reads_values = list(num_reads.values())
sort_cache_values = list(sort_cache.values())
scache_err_values = list(scache_err.values())
parse_cache_values = list(parse_cache.values())
pcache_err_values = list(pcache_err.values())

# Fig: cache miss values for Phase 1 and Phase 2 to validate our model
fig, axs = plt.subplots(1, 2, figsize=(6, 3))
for i, ax in enumerate(axs.flat):
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)

sns.lineplot(x=num_reads_values, y=parse_cache_values, marker='o', label='PAPI', legend=False, ax=axs[0], color=colors[1])
sns.lineplot(x=num_reads_values, y=estparsecache, marker='s', label='Model', legend=False, ax=axs[0], color=colors[0])
axs[0].errorbar(num_reads_values, parse_cache_values, xerr=pcache_err_values, fmt='', color=colors[1], ecolor=colors[1], capsize=5)
axs[0].set_title("Phase 1")
axs[0].set_xlabel("Number of Reads")
axs[0].set_ylabel("LL Cache Misses")
axs[0].tick_params(axis='both', which='major', labelsize=8)


sns.lineplot(x=num_reads.values(), y=sort_cache_values, marker='o', label='PAPI', ax=axs[1], color=colors[1])
sns.lineplot(x=num_reads.values(), y=estsortcache, marker='s', label='Model', ax=axs[1], color=colors[0])
axs[1].errorbar(num_reads_values, sort_cache_values, xerr=scache_err_values, color=colors[1], fmt='', ecolor=colors[1], capsize=5)
axs[1].legend()
axs[1].set_title("Phase 2")
axs[1].set_xlabel("Number of Reads")
axs[1].tick_params(axis='both', which='major', labelsize=8)

plt.tight_layout()

# Save the plot as an SVG file
plt.savefig('figures/hwcounts.svg', format='svg')
# plt.show()

# Fig: Time taken in Phase 1 and Phase 2 as predicted by our model and measured
fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharex=True, sharey=True)

for i, ax in enumerate(axs.flat):
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)

sns.lineplot(x=num_reads.values(), y=p1time_sum, marker='o', label='Model (Sum)', ax=axs[0], color=colors[0])
sns.lineplot(x=num_reads.values(), y=p1time_max, marker='*', label='Model (Max)', ax=axs[0], color=colors[1])
sns.lineplot(x=num_reads.values(), y=p1_time.values(), marker='s', label='Experiment', ax=axs[0], color=colors[2])
axs[0].set_title("Phase 1")
axs[0].set_xlabel("Number of Reads")
axs[0].set_ylabel("Time (seconds)")
axs[0].tick_params(axis='both', which='major', labelsize=8)

sns.lineplot(x=num_reads.values(), y=p2time, marker='o', label='Model', ax=axs[1], color=colors[0])
sns.lineplot(x=num_reads.values(), y=p2_time.values(), marker='s', label='Experiment', ax=axs[1], color=colors[2])
axs[1].legend()
axs[1].set_title("Phase 2")
axs[1].set_xlabel("Number of Reads")
axs[1].tick_params(axis='both', which='major', labelsize=8)

plt.tight_layout()

plt.savefig('figures/modeltime.svg', format='svg')
# plt.show()