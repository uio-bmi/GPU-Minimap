#ifndef MM_CUDA_H
#define MM_CUDA_H
#endif

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include "kseq.cuh"
#include <time.h>

/* Thrust */
#include <thrust/remove.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

/* CUB */
#include <cub/cub.cuh>
#include <cub/util_allocator.cuh>

#define NSTREAMS 32

typedef uint64_t minimizer_t;

typedef struct {
  minimizer_t m;
  int length;
} A_t;

/* device pointer, used for dna string on device */
extern int n_dna_memory;
extern char * d_dna;


/* Holds kmers. Approx 4 times the length of the sequence */
extern minimizer_t * d_kmers;
extern int n_kmer_memory;

extern minimizer_t * d_storage_mem;
extern int n_storage_memory;

extern minimizer_t * d_index_mem;
extern int n_index_memory;

extern minimizer_t * d_window_mem;
extern minimizer_t * d_window_ptr;
extern minimizer_t * d_query_mem;
extern minimizer_t * d_A;

extern int n_query_mem;
extern int n_A_mem;
extern int n_window_memory;
extern int global_window_end;

extern minimizer_t * overlap_set;

extern minimizer_t * cub_temp_mem[NSTREAMS];

extern size_t n_cub_mem;

#define CUDA_DEVICE_MAX_CONNECTIONS NSTREAMS

/* For lookup-table */
extern unsigned char * d_table;

/* All processing memory */
extern char * d_proc_mem;
extern int n_proc_sum;

#define WARP_SIZE 32
#define THREADS_PER_BLOCK 1024 // 768 on < 3.x devices.
#define WARPS_PER_BLOCK THREADS_PER_BLOCK / WARP_SIZE
#define TABLE_SIZE 256

/* Will organize all files when I'm done :p */

typedef struct {

    int k;    // K-mer
    int w;    // window size
    int r;    // (r) How large difference in diagonals allowed (500 default)
    float m;  // minimizeroverlap % (default 0.50 (50%)); What does this mean ?: "Minimap merges two windows if -m=50% of minimizer matches overlap."
    int c;    // minimum matches (chain length)
    int L;    //matching length
    int f;    // N/A
    int g;    //Break the chain whenever there is a gap longer than -g=10000.
    int n;    // N/A

} options_t;

typedef struct {

  /* Value */
  uint32_t min_value;
  uint32_t index;

} value_t;

typedef struct entry_t entry_t;

struct entry_t {

  /* Key */
  uint32_t key;

  /* Value */
  value_t val;

  /* Count */
  uint64_t count;


};

/* Section size */
#define OVERLAP_BLK (1LL << 27)

typedef struct {

  minimizer_t * set;
  int offset[NSTREAMS];
  int limit[NSTREAMS];

} overlap_set_t;

/* info about chromosome/read */
typedef struct info {

    char * name;
    int length;

} info_t;

/* Predicates */
struct is_ambiguous
{

    __host__ __device__ __forceinline__
    bool operator() (const uint64_t &m) {
      return(m == ((uint64_t)-1));
    }
};

struct AmbiguousStrand
{

    CUB_RUNTIME_FUNCTION __host__ __device__ __forceinline__
    bool operator() (const minimizer_t &m) {
      return (m != ((minimizer_t)-1));
    }
};

struct LessThan
{
  int compare;

  CUB_RUNTIME_FUNCTION __forceinline__
  LessThan(int compare) : compare(compare) {}

    CUB_RUNTIME_FUNCTION __forceinline__ __host__ __device__
    bool operator()(const int &a) const {
        return (a < compare);
    }
};

struct index_search {

  __host__ __device__
  bool operator() (const minimizer_t &m1, const minimizer_t &m2) {

    return ((m1 >> 32) < (m2 >> 32));
  }
};

void init_nres();

__global__
void encode_string(char * seq, unsigned char * table, minimizer_t * res, const int k, const int w, int rid, int SHIFT, int BIT_MASK, int length, int offset);

__global__
void find_kmers(char * seq, minimizer_t * res, const int k, const int w, int rid, int SHIFT, int BIT_MASK, int length, int offset);

__device__
inline minimizer_t cmp_min(minimizer_t m1, minimizer_t m2);

__global__ /* Will need to do some experimentation on whether all this shared memory fuzz is worth it in the end.. */
void reduce_to_minimizers(minimizer_t * minimizers, minimizer_t * results, int w, int length, int offset);

__global__
void find_overlaps(minimizer_t * A, minimizer_t * queries, minimizer_t * windows, minimizer_t * indices, unsigned int* lookup_indices, int index_len, int n, int * nres);

__global__
void compact_index(minimizer_t * index, minimizer_t * mm, int * nres, int n);

__host__ __device__
inline int get_index(minimizer_t m);

__host__ __device__
inline int get_strand(minimizer_t m);

__host__ __device__
inline int get_minimizer_value(minimizer_t m);

cudaError_t checkCuda(int line, cudaError_t result);

int calc_grid(int len, int n);

#define DEV_MEMORY (1 << 32) + (1 << 30)
#define DEV_BLOCK_SIZE (1 << 26)
#define BIAS (1 << 28)
#define DEV_BLOCK_NUM DEV_MEMORY / DEV_BLOCK_SIZE /* 64 (24)  */
#define MAX_MINIMIZER_OCC 8 // Maximum number of minimizer hits.. tunable

void setup_proc_memory(cudaStream_t * streams);
void free_proc_memory(cudaStream_t * streams);
void setup_index_memory();
void free_index_memory();
