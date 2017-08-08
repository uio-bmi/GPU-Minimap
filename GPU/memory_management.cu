#ifndef MEM_MGMT
#define MEM_MGMT
#endif

#include <time.h>

#include "minimap_cuda.cuh"

/* Alloc globally on device */
/* From lh3/minimap/sketch.c */
/* Rework */
unsigned char h_seq_table[256] = {
  0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/* CUB allocator (for temp) */
cub::CachingDeviceAllocator temp_allocator(true);

void setup_proc_memory(cudaStream_t * streams) {

    n_proc_sum = (DEV_BLOCK_SIZE * 9) + TABLE_SIZE; // 1 byte for char, 4 byte for kmers.

    n_dna_memory = DEV_BLOCK_SIZE;
    cudaMalloc((void**)&d_dna, n_dna_memory);

    n_kmer_memory = DEV_BLOCK_SIZE << 3;
    cudaMalloc((void**)&d_kmers, n_kmer_memory);

    cudaMalloc((void**)&d_table, 256);
    cudaMemcpy(d_table, h_seq_table, 256*sizeof(unsigned char), cudaMemcpyHostToDevice);
    /* d_table should end at the n_proc_sum byte */


    if ((overlap_set = (minimizer_t*)malloc(OVERLAP_BLK<<5)) < 0) {
        fprintf(stderr, "Error allocating overlap set\n");
    }

    checkCuda(__LINE__, cudaHostRegister((void*)overlap_set, OVERLAP_BLK<<5, 0));

    n_cub_mem = ((1<<16)<<4);

    for (int i = 0; i < NSTREAMS; i++) {
        checkCuda(__LINE__, temp_allocator.DeviceAllocate((void**)&cub_temp_mem[i], n_cub_mem, streams[i]));

    }

}

void setup_index_memory() {

     /* spare memory used for index, */
     n_storage_memory = DEV_BLOCK_SIZE << 2; //268435456

     /* Index table memory */
     n_index_memory = n_storage_memory<<1;
     cudaMalloc((void**)&d_index_mem, n_index_memory);//53687092

     /* Minimizer table */
     n_window_memory = n_storage_memory << 2;
     cudaMalloc((void**)&d_window_mem, n_window_memory);//107374184

     /* Readjust to one big malloc ^v */

     /* Overlap query mem */
     n_query_mem = DEV_BLOCK_SIZE<<3;
     cudaMalloc((void**)&d_query_mem, n_query_mem);


     /* A pointer to within windows */
     //d_window_ptr = &d_storage_mem[n_index_memory];
     d_window_ptr = d_window_mem;


}

void free_proc_memory(cudaStream_t * streams) {

    cudaFree(d_dna);
    cudaFree(d_kmers);
    cudaFree(d_table);
    cudaFree(d_query_mem);
    cudaFreeHost(overlap_set);


    for (int i = 0; i < NSTREAMS; i++) {
        temp_allocator.DeviceFree(cub_temp_mem[i]);
    }
}

void free_index_memory() {

    cudaFree(d_index_mem);
    cudaFree(d_window_mem);

}
