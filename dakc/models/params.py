# Phoenix parameters
m = 150
k = 31

# Phoenix parameters
blink = 12.5 * 10**9 # 12.5 * 10**9
bpcie = 64 * 10**9
A = 1 # 5.7 # Experimentally found

# CPU parameters
Z_cpu = 38 * 10**6 #/ 8  # 1 word = 8 bytes
L_cpu = 64 #/ 8          # 1 word = 8 bytes in Xeon CPUs
bmem_cpu = 47.16 * 10**9
cnode_cpu = 121.9 * 10**9 # 2.7 * 10**12

# GPU parameters
Z_gpu = 50 * 10**6 #/ 8
L_gpu = 128 #/ 8
bmem_gpu = 3 * 10**12
cnode_gpu = 25.6 * 10**12

# For isopower operation
tdp_cpu = 250 # For 2 x Xeon CPUs on Phoenix
tdp_gpu = 350
factor = tdp_gpu / tdp_cpu