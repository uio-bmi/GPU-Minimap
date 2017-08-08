#ifndef MM_CUDA
#define MM_CUDA
#endif

#include "minimap_cuda.cuh"

__device__
inline int invertible_hash(int x, int m) {

  x = (~x + (x << 21)) & m;
  x = x ^ x >> 24;
  x = (x + (x << 3) + (x << 8)) & m;
  x = x ^ x >> 14;
  x = (x + (x << 2) + (x << 4)) & m;
  x = x ^ x >> 28;
  x = (x + (x << 31)) & m;
  return x;

}

/* this kernel encodes string into 2-bit dna, calculates hashed kmers and stores in res */
__global__
void encode_string(char * seq, unsigned char * table, minimizer_t * res, const int k, const int w, int rid, int SHIFT, int BIT_MASK, int length, int offset)
{
  /* global index (into stream) */
  uint32_t g_idx = blockDim.x * blockIdx.x + threadIdx.x;

  /* aligned index (across all streams) */
  uint32_t a_idx = g_idx + offset;

  /* block index (index used for shared memory and within block) */
  uint32_t s_idx = threadIdx.x;

  /*
  * Shared
  */
  extern __shared__ unsigned char s_seq[];


  /* Return */
  if (g_idx >= length) return;

  /* Some files code N for */
  s_seq[s_idx] = (uint8_t)table[seq[g_idx]];

  /*First k threads handle overlaps (it is highly unlikely that k > # threads in warp (32) */

    if (s_idx < k - 1) {
      s_seq[s_idx + THREADS_PER_BLOCK] = (uint8_t)table[seq[g_idx + THREADS_PER_BLOCK]];
    }


  /* need seperate for last one */
  if (g_idx >= length - (k - 1)) return;

  __syncthreads();

  int i;
  int kmer[2] = {0,0};
  int c;
  uint32_t z;

  /* Make k-mer */
  #pragma unroll 15
  for (i = 0; i < k; i++) {

    c = s_seq[s_idx+i];

    /* skip any Ambiguous k-mers */
    if (c > 3) {
        res[a_idx] = (minimizer_t)-1;
        return;
    }

    kmer[0] = (kmer[0] << 2 | c) & BIT_MASK;           // forward k-mer
    kmer[1] = (kmer[1] >> 2) | (3UL^c) << SHIFT; // reverse k-mer


  }

 /* Symmetric kmers are eliminated */
if (kmer[0] == kmer[1]) {

    /* try to improve this */
    res[a_idx] = (minimizer_t)-1;
    return;
}

z = kmer[0] < kmer[1] ? 0 : 1; // strand

  uint32_t x = invertible_hash(kmer[z], BIT_MASK);

  /* TODO: set seq in here too */
  res[a_idx] = (((minimizer_t) x << 32) | (((rid << 29) | (a_idx)) << 1) | z);

}

/* get the minimizer value */
__device__ __host__
inline int get_minimizer_value(minimizer_t m) {

    return (m >> 32);
}

__device__ __host__
inline int get_index(minimizer_t m) {

    return ((m << 36)>>37);
}

__device__ __host__
inline int get_rid(minimizer_t m) {

    return ((m << 32) >> 61);
}
__device__ __host__
inline int get_strand(minimizer_t m) {

    return m & 1;
}

/* Might be better with cuda builtin*/
__device__
inline minimizer_t cmp_min(minimizer_t m1, minimizer_t m2) {

    return get_minimizer_value(m1) < get_minimizer_value(m2) ? m1 : m2;
}

__global__ /* Will need to do some experimentation on whether all this shared memory fuzz is worth it in the end.. */
void reduce_to_minimizers(minimizer_t * minimizers, minimizer_t * results, int w, int length, int offset) {

  /* this is the index into the global data */
  uint32_t aligned_idx = (blockIdx.x * offset) + threadIdx.x;

  /* block index (index used for shared memory) */
  uint32_t s_idx = threadIdx.x;


  /* Shared mem */
  extern __shared__ minimizer_t s_res[];


  /* Get from global */
  s_res[s_idx] = minimizers[aligned_idx];

  /* This part is fine, as we are assured that all memory was written to global in previous kernel. */
  if ((s_idx) < (w >> 1)) {

    s_res[s_idx + THREADS_PER_BLOCK] = minimizers[aligned_idx + THREADS_PER_BLOCK];

  }

  /*
   * If its the last block we must calculate the needed threads.
   */
  int s = (gridDim.x == blockIdx.x + 1) ? ((length % offset) - (w >> 1)) : THREADS_PER_BLOCK; // = len


  /* Make sure all threads are ready */
  __syncthreads();

  int r;
  int wc = w;

  minimizer_t temp;

  /* In the case of w=8, this will result in (8/2=4, 4/2=2, 2/2=1) -> 3 runs in this loop */
  /* In the case of w=5, this will result in (5/2=2, 3/2=1, 2/2=1) -> 3 runs in this loop aswell */
  #pragma unroll 4
  for (r = wc >> 1; r > 0; r >>= 1) {

    if (s_idx >= s) return;

    temp = cmp_min(s_res[s_idx], s_res[s_idx + r]); //res[p].m.x < res[g_idx + r].m.x ? p : g_idx + r;

    __syncthreads();    /* <-------------------- */
                                             /*| */
    s_res[s_idx] = temp;                     /*-------- need these to avoid race condition... ..*/
                                             /*| */
    __syncthreads();    /* <-------------------- */

    /*
    * readjust for odd numbers
    * For odd w, we must add +1. F.ex, 5 is not a factor of two, so we need (2*2+1)
    */
    r += wc & 1;
    wc = r;

    /* Reset thread threshold to be consistant with reduce offset */
    /* Important that this is done here and not before r update */
    s -= (r >> 1);

  }

  /* Write back results */
  results[aligned_idx] = s_res[s_idx];
}

/*
 * Copypasta (more or less) from https://devblogs.nvidia.com/parallelforall/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/
 * |
 * |
 * |
 * \/
 */

__device__
inline int lane_id(void) {
  return threadIdx.x % WARP_SIZE;
}

__device__
int warp_bcast(int v, int leader) {

  return __shfl(v, leader);
}

__device__
int atomicAggInc(int * ctr) {

  int mask = __ballot(1); // Active lanes

  int leader = __ffs(mask) - 1;

  int res;

  if (lane_id() == leader) {

    res = atomicAdd(ctr, __popc(mask));
  }

  res = warp_bcast(res, leader);

  return res + __popc(mask & ((1 << lane_id()) - 1));
}

__global__
void compact_index(minimizer_t * index, minimizer_t * mm, int * nres, int n) {

  /* Use shared ?*/

  /* Global thread index */
  int g_idx = blockDim.x * blockIdx.x + threadIdx.x;

  /* Out of range */
  if (g_idx >= n) return;

  /* If first */
  if (g_idx == 0) {
    index[atomicAggInc(nres)] = ((minimizer_t)get_minimizer_value(mm[g_idx]) << 32) | g_idx;
    return;
  }

  /* Add on predicate */
  if (get_minimizer_value(mm[g_idx - 1]) != get_minimizer_value(mm[g_idx])) {
    index[atomicAggInc(nres)] = ((minimizer_t)get_minimizer_value(mm[g_idx]) << 32) | g_idx;
  }

}

__global__ void find_overlaps(minimizer_t * A, minimizer_t * queries, minimizer_t * windows, minimizer_t * indices, unsigned int * lookup_indices, int index_len, int n, int * nres) {


  /* global index (into stream) */
  uint32_t g_idx = blockDim.x * blockIdx.x + threadIdx.x;

  /* block index (index used for shared memory and within block) */
  uint32_t t_idx = threadIdx.x;

  if (g_idx > n) {
    return;
  }

  minimizer_t q = queries[g_idx];

  unsigned int m = lookup_indices[g_idx];

  if (m >= index_len) {
    return;
  }


  int32_t pos = indices[m]<<32>>32;
  minimizer_t target;
  int bit;
  uint32_t diag;

  /* add > len too */
  if (pos < 0 ) {
    return;
  }

  #pragma unroll
  for (int i = 0; i < MAX_MINIMIZER_OCC; i++) {

    target = windows[pos+i];

    /* q, t don't match */
    if (target>>32 != q>>32) {
        return;
    }


    /* Same strand ? */
    /* Divergence on strands may cause some unwanted performance loss due to 2^2 combinations here.. */
    if (!((get_strand(q)) ^ get_strand(target))) {

      diag = (get_index(q) - get_index(target)) + BIAS;
      bit = 0x0;

    } else {

      diag = (get_index(q) + get_index(target)) + BIAS;
      bit = 0x1;
    }

     int ct = atomicAggInc(nres);

     uint64_t temp = ((uint64_t)get_rid(target));
     temp = temp << 1 | ((uint64_t)bit);
     temp = temp << 28 | ((uint64_t)diag);
     temp = temp << 28 | ((uint64_t)get_index(target));
     A[ct] = temp;

    }

}
