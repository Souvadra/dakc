import numpy as np
from math import log2, ceil

# log(N) base Z
def logZ(Z, N):
  return (log2(N)/log2(Z))

# kmernum, kmerbits and kmerbytes
def kmerbits(k):
  return 2**(ceil(log2(2 * k)))

def kmerbytes(k):
  return kmerbits(k) / 8

def kmernum(m, n, k): # n reads of m characters each
  return (n * (m - k + 1))

# N = Size of the intermediate array storing all the kmers per processor
def N(m, n, k, P): # also called all_kmer_size
  return kmernum(m, n, k) * kmerbytes(k) / P

def inputsize(m, n, P): # Output is in Bytes
  return (m * n) / P

# sorting complexity (hybrid sorting)
def sortcomp(x, k):
  # x is the number of kmers to sort 
  # each kmer needs kmerbytes(k) bytes
  # hence total bytes to sort = x * kmerbytes(k)
  radixsort = x * kmerbytes(k)
  return radixsort

# number of cache misses in the hybrid sorting
def sortcomm(x, Z, L, k, A2):
  # x -> amount of data in bytes 
  # radixsort will parse the data kmerbytes(k) times
  # each parsing causes (1 + (x / L)) cache misses
  # hence total cache misses = A2 * kmerbytes(k) * (1 + (x / L))
  return (A2 * kmerbytes(k) * (1 + (x / L)))

# Phase 1 computation:
def t1comp(m, n, P, cnode):
  return (inputsize(m, n, P) / cnode)

def t1intracachemiss(m, n, P, L, k):
  term1 = (1 + (inputsize(m, n, P) / L))
  term2 = (1 + (N(m, n, k, P) / L))
  assert(term1 > 0)
  assert(term2 > 0)
  return (term1 + term2)

def t1intra(m, n, P, L, bmem, k):
  term1 = t1intracachemiss(m, n, P, L, k)
  term2 = L / bmem
  return (term1 * term2)

def t1inter(m, n, k, P, blink):
  constant = 2
  # Constant is 2 because each processor is sending and receiving the same 
  # amount of data, and blink is bidirectional bandwidth
  return constant * (N(m, n, k, P) / blink)

# Phase 2 computation
def t2comp(m, n, k, P, cnode):
  datasize = kmernum(m, n, k) / P
  return (sortcomp(datasize, k) / cnode)

def t2intra(m, n, k, P, L, Z, bmem, A):
  datasize = N(m, n, k, P)
  term1 = L / bmem
  return (sortcomm(datasize, Z, L, k, A) * term1)

def cpu_time(m, n, k, P, L, Z, cnode, bmem, blink, A, type):
  t1_comp   = t1comp(m, n, P, cnode)
  t1_intra  = t1intra(m, n, P, L, bmem, k)
  t1_inter  = t1inter(m, n, k, P, blink)
  t2_comp   = t2comp(m, n, k, P, cnode)
  t2_intra  = t2intra(m, n, k, P, L, Z, bmem, A)

  t1_comm = 0
  t2_comm = 0
  if (type == "sum"):
      t1_comm += t1_intra + t1_inter
      t2_comm += t2_intra
  else:
      assert(type == "max")
      t1_comm += max(t1_intra, t1_inter)
      t2_comm += t2_intra

  t1 = max(t1_comp, t1_comm)
  t2 = max(t2_comp, t2_comm)

  return (t1 + t2)

def cpu_phase1time(m, n, k, P, L, Z, cnode, bmem, blink, A, type):
  t1_comp   = t1comp(m, n, P, cnode)
  t1_intra  = t1intra(m, n, P, L, bmem, k)
  t1_inter  = t1inter(m, n, k, P, blink)

  t1comm = 0
  if (type == "sum"):
      t1comm += t1_intra + t1_inter
  else:
      assert(type == "max")
      t1comm += max(t1_intra, t1_inter)

  t1 = max(t1_comp, t1comm)
  return t1

def cpu_phase2time(m, n, k, P, L, Z, cnode, bmem, blink, A):
  t2_comp   = t2comp(m, n, k, P, cnode)
  t2_intra  = t2intra(m, n, k, P, L, Z, bmem, A)
  t2 = max(t2_comp, t2_intra)
  return t2

def cpu_compute(m, n, k, P, L, Z, cnode):
  return (t1comp(m, n, P, cnode) + t2comp(m, n, k, P, cnode))

def cpu_intranode(m, n, k, P, L, Z, bmem, A):
  t1term = t1intra(m, n, P, L, bmem, k)
  t2term = t2intra(m, n, k, P, L, Z, bmem, A)
  return (t1term + t2term)

def cpu_internode(m, n, k, P, blink):
  return t1inter(m, n, k, P, blink)
