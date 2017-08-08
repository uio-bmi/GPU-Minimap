#ifndef MAIN
#define MAIN
#endif

#include "minimap_cuda.cuh"
#include <time.h>
#include <cuda_profiler_api.h>

KSEQ_INIT(int, read)

/* Alloc globally on device */
/* From lh3/minimap/sketch.c */
unsigned char hh_seq_table[256] = {
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


/* device pointer, used for dna string on device */
int n_dna_memory;
char * d_dna;


/* Holds kmers. Approx 4 times the length of the sequence */
minimizer_t * d_kmers;
int n_kmer_memory;

minimizer_t * d_storage_mem;
int n_storage_memory;

minimizer_t * d_index_mem;
int n_index_memory;

minimizer_t * d_window_mem;
minimizer_t * d_window_ptr;

minimizer_t * d_query_mem;
minimizer_t * d_A;

cudaStream_t streams[NSTREAMS];

minimizer_t * overlap_set;
minimizer_t * cub_temp_mem[NSTREAMS];

size_t n_cub_mem;

/*
 * Some restrictions to read lengths.
 * Haven't been thorougly tested for very short read lengths.
 * Should be fine, but exists to prevent unforeseen edge cases.
 * very short reads.
 */
#define MIN 2000
#define MAX 80000


int n_query_mem;
int n_A_mem;

int n_window_memory;
int global_window_end;

/* For lookup-table */
unsigned char * d_table;

/* All processing memory */
char * d_proc_mem;
int n_proc_sum;

/* Index size */
int h_nres;

int SHIFT;
int BIT_SIZE;
int BIT_MASK;

int calc_grid(int len, int n) {

  return (len / n) + ((len % n) != 0 ? 1 : 0);
}


/* copied from https://github.com/parallel-forall/code-samples/blob/master/series/cuda-cpp/overlap-data-transfers/async.cu */
/* For error checking */
cudaError_t checkCuda(int line , cudaError_t result)
{
  #ifdef DEBUG_CUDA

  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s, at line: %d\n", cudaGetErrorString(result), line);
    assert(result == cudaSuccess);
  }

  #endif

  return result;
}

int index_targets(char * path, int k, int w) {

  int k_overlap = k - 1;
  int w_overlap = w - 1;

  int stream_size;
  int rem;
  int n_grid;
  int last_grid;

  const int SHIFT = k_overlap << 1;
  const int BIT_SIZE = k << 1;
  const int BIT_MASK = ((int)(1 << BIT_SIZE)-1);

  const int BLOCK_OFFSET = (w >> 1) - w_overlap; /* Offset within block (always negative);*/
  const int RESULT_THREADS = THREADS_PER_BLOCK + BLOCK_OFFSET; /* Number of actual result yielding threads per block, given THREADS_PER_BLOCK number of threads.*/
  const int BLOCK_SIZE = THREADS_PER_BLOCK; // = 1024
  int INDEX_STREAMS = NSTREAMS;

  thrust::device_ptr<minimizer_t> dev_ptr_kmer_set = thrust::device_pointer_cast(d_kmers);
  thrust::device_ptr<minimizer_t> dev_ptr_window_begin = thrust::device_pointer_cast(d_window_mem);

  FILE * fp = fopen(path, "r");

  if (fp == 0) {

      fprintf(stderr, "Error with file pointer..\n");
      exit(0);
  }

  kseq_t * seq = kseq_init(fileno(fp));

  int c;
  int ct = 0;
  int len;
  int num_kmers, num_windows;

  /* Read sequences from file */
  while((c = kseq_read(seq)) >= 0) {

      len = seq->seq.l;
      if (len < 1000) continue;

      /* Waste using too many streams on small inputs */
      if (len < 50000000) {
          INDEX_STREAMS = 1;
      } else {
        INDEX_STREAMS =  (NSTREAMS < 4) ? NSTREAMS : 4;

      }

      /* Temporary bounds check for window set */

      /* Pin host memory */
      cudaHostRegister(seq->seq.s,len, 0);

      /* section is how much unique data processed per dev_block */
      int section = DEV_BLOCK_SIZE;
      int section_offset = section - k_overlap - w_overlap;

      int rounds = calc_grid(len, section_offset);

      /* The remaining part for last round */
      int round_rest = len % section_offset;

      for (int i = 0; i < rounds; i++) {

          /* Offset processed by each round */
          int memory_offset = i * section_offset;

          /* if last round (or first, if len < DEV_BLOCK_SIZE.. */
          if (i == rounds-1) {

              section = round_rest;

          }

          /* Number of k-mers section yields.. */
          num_kmers = section - k_overlap;


          /* We split each round into streams for better data transfer */
          stream_size = section / (INDEX_STREAMS); // try this for the time being.
          rem = section % (INDEX_STREAMS);  // ---- num_kmers % (INDEX_STREAMS)? ^
          n_grid = calc_grid(stream_size + k_overlap, THREADS_PER_BLOCK);
          last_grid = calc_grid(stream_size + rem, THREADS_PER_BLOCK);

          /* Unroll if stream count is constant */
          for (int j = 0; j < (INDEX_STREAMS); j++) {
              int stream_offset = j * stream_size;

              /* extra + k_overlap for overlap */
              if (j == (INDEX_STREAMS) -1) {

                /* For remainer (given a number not divisible with (INDEX_STREAMS)) */
                cudaMemcpyAsync(&d_dna[stream_offset], &(seq->seq.s[stream_offset + memory_offset]), stream_size + rem, cudaMemcpyHostToDevice, streams[j]);
                encode_string<<<last_grid, BLOCK_SIZE, BLOCK_SIZE + k_overlap, streams[j]>>>(&d_dna[stream_offset], d_table, d_kmers, k, w, ct, SHIFT, BIT_MASK, stream_size + rem, stream_offset);

              }
              else {

                cudaMemcpyAsync(&d_dna[stream_offset], &(seq->seq.s[stream_offset + memory_offset]), stream_size + k_overlap, cudaMemcpyHostToDevice, streams[j]);
                encode_string<<<n_grid, BLOCK_SIZE, BLOCK_SIZE + k_overlap, streams[j]>>>(&d_dna[stream_offset], d_table, d_kmers, k, w, ct, SHIFT, BIT_MASK, stream_size + k_overlap, stream_offset);

              }
          }


          /* We need this for our reduce method to work. */
          thrust::device_ptr<minimizer_t> end_pos_ambiguous = thrust::remove_if(thrust::device, dev_ptr_kmer_set, dev_ptr_kmer_set+num_kmers, is_ambiguous());


          /* Length after filtering out N's and ambiguous kmers */
          int filtered_length = end_pos_ambiguous - dev_ptr_kmer_set;

          /* Number of windows possible */
          num_windows = filtered_length - (w_overlap);

          /* Window offset already taken into account in result threads */
          const int WINDOW_GRID_SIZE = calc_grid(filtered_length, RESULT_THREADS);
                                                                              /* +1?*/
          reduce_to_minimizers<<<WINDOW_GRID_SIZE, BLOCK_SIZE, (BLOCK_SIZE + (w >> 1)) * sizeof(minimizer_t)>>>(d_kmers, d_window_ptr, w, filtered_length, RESULT_THREADS);

          /* At this point, all minimizers should be found and we need only retrieve the unique */
          thrust::device_ptr<minimizer_t> dev_ptr_window_set = thrust::device_pointer_cast(d_window_ptr);


          /* Move out? */
          thrust::device_ptr<minimizer_t> end_pos = i == 0 ? thrust::unique(thrust::device, dev_ptr_window_set, dev_ptr_window_set+num_windows) : thrust::unique(thrust::device, dev_ptr_window_set-1, dev_ptr_window_set+num_windows);

          /* where the set ends (for sorting etc.) */
          global_window_end = end_pos - dev_ptr_window_begin - 1; /* -1 because end_pos is where the undefined part starts */

          /* Set pointer to end */
          d_window_ptr = thrust::raw_pointer_cast(end_pos);
          //global_window_end = num_windows;
          /* Also need ot check that pointer doesnt exceed limits etc. */
      }

      /* Unpin */
      cudaHostUnregister(seq->seq.s);
  }

  /* Sorting phase */

  thrust::sort(thrust::device, dev_ptr_window_begin, dev_ptr_window_begin + global_window_end);
  thrust::device_ptr<minimizer_t> dev_ptr_index = thrust::device_pointer_cast(d_index_mem);

  int compact_grid = calc_grid(global_window_end, THREADS_PER_BLOCK);
  int * d_nres;
  h_nres = 0;

  /* Can probably solve this better with __device__ variables .. */
  cudaMalloc((void**)&d_nres, sizeof(int));
  cudaMemcpy(d_nres, &h_nres, sizeof(int), cudaMemcpyHostToDevice);

  /* Filter minimizers into index table */
  compact_index<<<compact_grid, THREADS_PER_BLOCK>>>(d_index_mem, d_window_mem, d_nres, global_window_end);
  cudaMemcpy(&h_nres, d_nres, sizeof(int), cudaMemcpyDeviceToHost);

  thrust::sort(thrust::device, dev_ptr_index, dev_ptr_index + h_nres);

  /*
   * At this point our table is done and we can create our lookup table. Binary search vs. hashing (cuckoo) will be
   * discussed in thesis. BS was chosen due to simplicity and performance (src).
   */

   /* There are several solutions to implementing this, but right now I'm thinking
    * that copying unique (first) entries into the index table and then sorting it
    * is the best, and by far the easiest.
    *
    * thrust copyif is probably the best bet, as it is optimized well for this task
    * possibly use https://devblogs.nvidia.com/parallelforall/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/ too.
    * as it shows great promise on our working architecture, as we don't necessarily require our set to be ordered the same
    * way during the insertion (as we can sort the set after). possibly use lowerbound thrust function too.
    *
    * This leaves us with some options;
    *
    * 1. Cuckoo hashing based on http://idav.ucdavis.edu/~dfalcant//downloads/dissertation.pdf
    *     - This is probably the most advanced solution, but this is not necessarily justified by
    *       it's performance in my opinion, as my time is limited.
    *
    * 2. Sorted lookup table using binary search.
    *     - This is pretty simple, we can make use of the optimized filtering mentioned earlier,
    *       and most likely using thrust lower bound algorithm (binary search).
    *
    * 3. Binary search directly into window set (of minimizers).
    *     - In this solution, we don't create a seperate table, but use binary search
    *       directly into the window set. This is by far the easiest to implement,
    *       but seems like the least efficient in terms of performance. On the other
    *       hand, we reduce the overhead of constructing a hash/lookup-table, thus
    *       thus limiting the amount of time and memory used for this.
    *
    *
    * Generally, for all binary search implementations, the performance gain
    * for binary search is improved greatly by having ordered queries, which
    * leads to the question of whether it is a good idea to sort them before
    * lookup.
    */

  fclose(fp);

  return 0;

}

int map_queries(char * path, int k, int w) {

  cudaEvent_t start_event, end_event;

  const int D_KMER_OFFSET = (n_kmer_memory / NSTREAMS) >>3;
  const int D_QUERY_OFFSET = (n_query_mem / NSTREAMS) >>3;

  int k_overlap = k - 1;
  int w_overlap = w - 1;

  const int SHIFT = k_overlap << 1;
  const int BIT_SIZE = k << 1;
  const int BIT_MASK = ((int)(1 << BIT_SIZE)-1);

  const int BLOCK_OFFSET = (w >> 1) - w_overlap; /* Offset within block (always negative);*/
  const int RESULT_THREADS = THREADS_PER_BLOCK + BLOCK_OFFSET; /* Number of actual result yielding threads per block, given THREADS_PER_BLOCK number of threads.*/
  const int BLOCK_SIZE = THREADS_PER_BLOCK; // = 1024

  checkCuda(__LINE__, cudaEventCreate(&start_event));
  checkCuda(__LINE__, cudaEventCreate(&end_event));

  FILE * fp = fopen(path, "r");

  if (fp == 0) {

      fprintf(stderr, "Error with file pointer..\n");
      exit(0);
  }

  kseq_t * seq = kseq_init(fileno(fp));

  int c = 0;

  const static unsigned max_length = (1 << 30); // tunable

  char * seqs;
  int * seq_lengths;

  cudaMallocHost((void**)&seq_lengths, (1<<15)*sizeof(int));
  seqs = (char*)malloc(max_length);
  std::vector<int> start_pos;

  int seq_lengths_sum = 0;
  int read_counter = 0;


  /* Get reads */
  while (kseq_read(seq) >= 0) {

      /* Can tweak this, but defined by macros up top. */
      if (seq->seq.l < MIN || seq->seq.l > MAX) {
          continue;
      }

      start_pos.push_back(seq_lengths_sum);

      memmove(&seqs[seq_lengths_sum], seq->seq.s, seq->seq.l);
      seq_lengths_sum+=seq->seq.l;

      seq_lengths[read_counter++] = seq->seq.l;
      if (seq_lengths_sum >= max_length) {
          break;
      }


  }

  char * d_reads;
  checkCuda(__LINE__, cudaMalloc((void**)&d_reads, seq_lengths_sum));
  checkCuda(__LINE__, cudaHostRegister((void*)&seqs[0], seq_lengths_sum, 0));

  int filtered_length[NSTREAMS];

  int set_length[NSTREAMS];

  int transfer_offset = seq_lengths_sum / NSTREAMS;//seqs.size() / NSTREAMS;
  int transfer_rem = seq_lengths_sum % NSTREAMS;//seqs.size() % NSTREAMS;

  /* Could have overlapped some kernels here I guess.. */
  int nreads = read_counter;

  for (int i = 0; i < NSTREAMS; i++) {

      if (i < NSTREAMS-1) {
        checkCuda(__LINE__, cudaMemcpyAsync(&d_reads[i*transfer_offset], &seqs[i*transfer_offset], transfer_offset, cudaMemcpyHostToDevice, streams[i]));
        continue;
      }

      checkCuda(__LINE__, cudaMemcpyAsync(&d_reads[i*transfer_offset], &seqs[i*transfer_offset], transfer_offset+transfer_rem, cudaMemcpyHostToDevice, streams[i]));
  }

  int rounds = calc_grid(read_counter, NSTREAMS);

  int num_kmers;
  int num_windows;

  const int OVERLAP_BLK_OFFSET = OVERLAP_BLK>>3;

  int overlap_limit[NSTREAMS];
  int overlap_offset[NSTREAMS];

  /* For binary search results */
  thrust::device_vector<unsigned int> query_outputs[NSTREAMS];

  /* Init some stuff.. */
  for (int i = 0; i < NSTREAMS; i++) {

      query_outputs[i].resize(1<<16);
      overlap_offset[i] = 0;
      overlap_limit[i] = OVERLAP_BLK_OFFSET;
  }

  /* Lengths from cub outputs */
  int * n_cub_output;
  checkCuda(__LINE__, cudaMalloc((void**)&n_cub_output, NSTREAMS*sizeof(int)));

  /* transfer_counters keeps track of elements ready for transfer */
  int transfer_counters[NSTREAMS];

  /* overlap_counters keeps track of currently processed elements */
  int overlap_counters[NSTREAMS];

  checkCuda(__LINE__, cudaEventRecord(end_event, 0));
  checkCuda(__LINE__, cudaEventSynchronize(end_event));

  /* Init counters to zero */
  for (int i = 0; i < NSTREAMS; i++) {
      transfer_counters[i] = 0;
      overlap_counters[i] = 0;
  }

  /* Check done by cub::unique::if etc. */
  AmbiguousStrand check;

  for (int i = 0; i < rounds; i++) {

    for (int j = 0; j < NSTREAMS && j < rounds; j++) { // double check this TODO

      /* index into pointers/lengths etc. */
      int l_idx = j*rounds + i;

      if (l_idx >= nreads) {
          break;
      }

      /* Some indices */
      int len = seq_lengths[l_idx];
      int s_p = start_pos[l_idx];

      int grid_size = calc_grid(len, BLOCK_SIZE);

      /* offset in d_kmers or as offset parameter */
      encode_string<<<grid_size, BLOCK_SIZE, BLOCK_SIZE + k_overlap, streams[j]>>>(&d_reads[s_p], d_table, &d_kmers[j*D_KMER_OFFSET], k, w, 0, SHIFT, BIT_MASK, len, 0);


      /*
       * At this point we can either conclude the loop, waiting for all loops
       * to complete before proceeding, or we can continue asynchronous stream
       * operations as usual as this will not affect other the CPU/other streams
       */
       num_kmers = len - k_overlap;

       cub::DeviceSelect::If((void*)cub_temp_mem[j], n_cub_mem, &d_kmers[j*D_KMER_OFFSET], &d_kmers[j*D_KMER_OFFSET+(D_KMER_OFFSET>>1)], &n_cub_output[j], num_kmers, check, streams[j]);

     }

     /* copy filtered lengths for CUB call */
     checkCuda(__LINE__, cudaMemcpy(filtered_length, n_cub_output, sizeof(int)*NSTREAMS, cudaMemcpyDeviceToHost));

     for (int j = 0; j < NSTREAMS && j < rounds; j++) { // double check this TODO

       /* index into pointers/lengths etc. */
       int l_idx = j*rounds + i;

       if (l_idx >= nreads) {
           // printf("Done Overlapping ...\n");
           break;
       }

       /* Use filtered_length here, but could easily be done inside kernels */
       num_windows = filtered_length[j] - (w_overlap);
       const int WINDOW_GRID_SIZE = calc_grid(filtered_length[j], RESULT_THREADS);
       reduce_to_minimizers<<<WINDOW_GRID_SIZE, BLOCK_SIZE, (BLOCK_SIZE + (w >> 1)) * sizeof(minimizer_t), streams[j]>>>(&d_kmers[j*D_KMER_OFFSET+(D_KMER_OFFSET>>1)], &d_kmers[j*D_KMER_OFFSET], w, filtered_length[j], RESULT_THREADS);

       /* At this point, all minimizers should be found and we need only retrieve the unique */

       cub::DeviceSelect::Unique((void*)cub_temp_mem[j], n_cub_mem, &d_kmers[j*D_KMER_OFFSET], (minimizer_t*)&d_kmers[j*D_KMER_OFFSET+(D_KMER_OFFSET>>1)], &n_cub_output[j], num_windows, streams[j]);
     }

     checkCuda(__LINE__, cudaMemcpy(set_length, n_cub_output, sizeof(int)*NSTREAMS, cudaMemcpyDeviceToHost));

     for (int j = 0; j < NSTREAMS && j < rounds; j++) { // double check this TODO

           /* index into pointers/lengths etc. */
           int l_idx = j*rounds + i;

           if (l_idx >= nreads) break;

           minimizer_t * queries = (minimizer_t*)&d_kmers[j*D_KMER_OFFSET+(D_KMER_OFFSET>>1)];

           thrust::lower_bound(thrust::cuda::par.on(streams[j]), thrust::device_ptr<minimizer_t>(d_index_mem), thrust::device_ptr<minimizer_t>(d_index_mem+h_nres), thrust::device_ptr<minimizer_t>(queries), thrust::device_ptr<minimizer_t>(queries+set_length[j]), query_outputs[j].begin(), index_search());

          /* Locate overlaps */
           checkCuda(__LINE__, cudaMemsetAsync(&n_cub_output[j], 0, sizeof(int), streams[j]));
           find_overlaps<<<calc_grid(set_length[j], BLOCK_SIZE), BLOCK_SIZE, 0, streams[j]>>>((minimizer_t*)&d_kmers[j*D_KMER_OFFSET], queries, d_window_mem, d_index_mem, thrust::raw_pointer_cast(&((query_outputs[j])[0])), h_nres, set_length[j], &n_cub_output[j]);

           /* Set start pos of read l_idx */
           start_pos[l_idx] = overlap_counters[j];

           /* Reading in size of read result */
           checkCuda(__LINE__, cudaMemcpyAsync(&seq_lengths[l_idx], &n_cub_output[j], sizeof(int), cudaMemcpyDeviceToHost, streams[j]));
     }

     /* Synchronise transfers */
     checkCuda(__LINE__, cudaEventRecord(end_event));
     checkCuda(__LINE__, cudaEventSynchronize(end_event));

     for (int j = 0; j < NSTREAMS && j < rounds; j++) { // double check this TODO

       /* index into pointers/lengths etc. */
       int l_idx = j*rounds + i;

       if (l_idx >= nreads) {
           break;
       }

       /* Sort here? or block sort at end of one loop. - must sync + more work per dispatch. */


       /* Room for more data in current temp memory block ? */
       if ((D_QUERY_OFFSET - transfer_counters[j]) < seq_lengths[l_idx]) {

            /* check if remaning is less than required */
            if (overlap_limit[j] - overlap_offset[j] < transfer_counters[j]) {
                printf("Need reallocation of %d's overlap set\n", j);
                exit(0);

            }

            /* Cross-check referencing here TODO */
            checkCuda(__LINE__, cudaMemcpyAsync(&(overlap_set[j*OVERLAP_BLK_OFFSET + overlap_offset[j]]), &d_query_mem[j*D_QUERY_OFFSET], sizeof(minimizer_t)*transfer_counters[j], cudaMemcpyDeviceToHost, streams[j]));
            overlap_offset[j] += transfer_counters[j];
            transfer_counters[j] = 0;
        }
    }

    for (int j = 0; j < NSTREAMS && j < rounds; j++) {

        /* index into pointers/lengths etc. */
        int l_idx = j*rounds + i;

        if (l_idx >= nreads) break;

             cub::DeviceRadixSort::SortKeys((void*)cub_temp_mem[j], n_cub_mem, (uint64_t*)&d_kmers[j*D_KMER_OFFSET], (uint64_t*)&d_query_mem[j*D_QUERY_OFFSET+transfer_counters[j]], seq_lengths[l_idx], 0, sizeof(uint64_t)*8, streams[j]);

             /* Increment transfer counter */
             transfer_counters[j] += seq_lengths[l_idx];

             /* Increment overlap counter */
             overlap_counters[j] += seq_lengths[l_idx];


          }

    }

/* Copy any (remaining) results over to host */
for (int i = 0; i < NSTREAMS; i++) {

    /* If any remaining */
    if (transfer_counters[i] > 0) {

        /* check if remaning is less than required */
        if (overlap_limit[i] - overlap_offset[i] < transfer_counters[i]) {

            fprintf(stderr, "Need reallocation of %d's overlap set\n", i);
            exit(0);

        }

        checkCuda(__LINE__, cudaMemcpyAsync(&overlap_set[(i*OVERLAP_BLK_OFFSET) + overlap_offset[i]], &d_query_mem[i*D_QUERY_OFFSET], sizeof(minimizer_t)*transfer_counters[i], cudaMemcpyDeviceToHost, streams[i]));

        overlap_offset[i] += transfer_counters[i];
        transfer_counters[i] = 0;
    }

}

  /* End event */
  checkCuda(__LINE__, cudaEventRecord(end_event, 0));
  checkCuda(__LINE__, cudaEventSynchronize(end_event));

  /* For checking number of overlaps */
  // int sum = 0;
  // for (int i = 0; i < NSTREAMS; i++) {
  //     sum += overlap_offset[i];
  // }
  // printf("num:%d\n", sum);

  /*Free up stuff */
  cudaFree(n_cub_output);
  cudaFree(d_reads);
  free(seqs);

  return 0;
}


/* argv 1 -> Ref, argv 2 -> Query, argv 3 -> k, argv 4 -> w */
int main(int args, char ** argv) {

  cudaEvent_t start, end;
  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

  const int k = atoi(argv[3]);
  const int w = atoi(argv[4]);
  int k_overlap = k - 1;

  SHIFT = k_overlap << 1;
  BIT_SIZE = k << 1;
  BIT_MASK = ((int)(1 << BIT_SIZE)-1);

  /* Init streams */
  for (int i = 0; i < NSTREAMS; i++)
      checkCuda(__LINE__,  cudaStreamCreate(&streams[i]));

  /* Init memory */
  setup_proc_memory(streams);
  setup_index_memory();

  /* Index */
  index_targets(argv[1], k, w);

  /* map_queries */
  map_queries(argv[2], k, w);

  /* Free up */
  free_proc_memory(streams);
  free_index_memory();

  for (int i = 0; i < NSTREAMS; i++)
    checkCuda(__LINE__, cudaStreamDestroy(streams[i]));


  /* Free memory */
  cudaEventDestroy(start);
  cudaEventDestroy(end);

}
