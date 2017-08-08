#ifndef MAP_CU
#define MAP_CU
#endif

#include "minimap_cuda.cuh"



inline uint32_t get_A_sequence(minimizer_t i) {

  return (uint32_t)((i >> 61));
}

inline uint32_t get_A_strand(minimizer_t i) {

  return (uint32_t)((i << 7) >> 63);
}

inline int32_t get_A_diagonal(minimizer_t i) {

  return (int32_t)((i << 8) >> 32);
}

inline uint32_t get_A_index(minimizer_t i) {

  return (uint32_t)(i<<36>>36);
}


/* Return the query index back from diagonal */
inline uint32_t convert_index(minimizer_t a) {

  if (get_A_strand(a) == 0) {

    return get_A_diagonal(a) + get_A_index(a);

  } else {

    return get_A_diagonal(a) - get_A_index(a);
  }
}

/* Used with qsort in map() */
inline int compare_query_index(const void * a, const void * b) {

		return convert_index(*(const minimizer_t*)a) - convert_index(*(const minimizer_t*)b);
}

inline void print_chain(minimizer_t * A, uint32_t * lis/*,kseq_t * q */, int b, int e) {
	int sequence = get_A_sequence(A[b]);
  //print overlap region
}

int map(unsigned long * A, int A_counter) {

  int b = 0;
  int e;

  int ct;


	int diff = 0;
	uint32_t * lis = (uint32_t*)malloc(1);
	uint32_t * p = (uint32_t*)malloc(1);

  /* Clusters within radius */
  for (e = 1; e < A_counter; e++) {
    if ((e == A_counter - 1)
    || (get_A_sequence(A[e + 1]) != get_A_sequence(A[e]))
    || (get_A_strand(A[e + 1]) != get_A_strand(A[e]))
    || (abs(get_A_diagonal(A[e + 1]) - get_A_diagonal(A[e]))) >= opts.r)
    {


      diff = e - b + 1; // +1 for 0-index
      /* Sort this interval on query index ...
       * We must sort because there may be hits in the cluster which have a lower diagonal, but higher q_index.
       * ie.
       *   _____________________________________________________________________
       *	|																						\
       *  |																						 \
       *  |																							\
       *  |																							 \
       *  |																							    \ <---- lower diag, higher q_index
       *
       *
       */




      /*
      *
      * C <-- The maximal colinear subset of A[b..e] (Longest Increasing Sequence Problem).
      * Print the left- and right-most query/target positions in C
      * b <-- e + 1
      *
      */

      /* We check this here aswell to avoid uneccessary sorting. */
      if (diff >= opts.c) {

        /* Lets try with qsort first. */
        qsort(&A[b], diff, sizeof(minimizer_t), compare_query_index);

        /* This should be where we find the longest increasing sequence. */

        if ((lis = (uint32_t*)realloc((void**)lis, sizeof(uint32_t)*diff)) == NULL) {
            fprintf(stderr, "Error allocating lis in minimap.c\n");
            exit(0);
        }

        if ((p = (uint32_t*)realloc((void**)p, sizeof(uint32_t)*diff)) == NULL) {
            fprintf(stderr, "Error allocating p in map.\n");
            exit(0);
        }

        int rev = get_A_strand(A[b]);

        int length_lis = longest_increasing_subsequence(&A[b], lis, p, rev, diff);
        int t;

        /* Proceed to break LIS into regions based on opts.g (gap) in query index. */

        /* Print regions if they satisfy >= c/L criteria */

        int inner_e = 1;
        int inner_b = 0;

        uint32_t val;

        if (length_lis >= opts.c) {

          while (inner_e < length_lis-1) {

              val = convert_index(A[b + lis[inner_e]]);

              /* Break chain if large enough gap */
              if (abs((int32_t)(convert_index(A[b + lis[inner_e + 1]]) - val)) >= opts.g) {

                  int temp = (inner_e - inner_b + 1);

                  /* If chain contains atleast -c=4 minimizers, and total length is atleast -L=40, then output chain. */
                  if (A[b + lis[inner_e]] - A[b + lis[inner_b]] + opts.k >=  opts.L && temp >= opts.c) { // -L fix parameter
                      
                      // print_chain(&A[b], lis, q, inner_b, inner_e);
                  }

                  inner_b = inner_e + 1;
              }

              inner_e++;
          }

          int temp = (inner_e - inner_b + 1);

          /*
           * Prints rest of chain if it satisfies the criterias
           */
          if (A[b + lis[inner_e]] - A[b + lis[inner_b]] + opts.k >= opts.L && temp >= opts.c) { // -L fix parameter

              // print_chain(&A[b], lis, q, inner_b, inner_e);
          }
        }

      }

      b = e + 1;

    }
  }

  free(lis);
  free(p);

  return 0;
}
