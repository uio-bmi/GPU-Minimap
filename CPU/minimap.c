#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include "minimap.h"

#define MIN(x, y) ((x < y) ? x : y)
#define MINIMIZER_BLOCK 1000000 // This should really scale to file sizes, but assuming files are very large.
#define BIAS (1 << 30)

/* From lh3/minimap/sketch.c */
unsigned char seq_nt4_table[256] = {
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

uint64_t BIT_SIZE;
uint64_t BIT_MASK;

/* Hash table */
hash_table_t * ht;

int overlap_count;

static options_t opts;


/* minimizers */
minimizer_set_t minimizer_target_set;
minimizer_set_t minimizer_query_set;

/* kseqs*/
kseq_t ** sequences;
int * seq_l;
kseq_t * map_sequence;

int overlap_count;
int main(int argc, char ** argv) {

	overlap_count = 0;

	opts = (options_t)   {.k = atoi(argv[3]),		// k-mer
												.w = atoi(argv[4]), 		// Window
												.r = 500,	// Radius
												.m = 0.5,		// Minimizer overlap
											  .c = 4,			// Chain length
											 	.L = 40,		// Window length //TODO WRONG
											  .f = 0.1, 	// N/A
											  .g = 10000, 	// Gap length
												.n = argc-4,
											};

	/* Index targets */
  int c = read_targets(argc, argv);


  if (c < 0) {
    fprintf(stderr, "Error handling file input\n");
    exit(0);
  }

  hash_table_t * ht = mm_index();

  int cr;

	/* Start mapping reads/queries */
  if ((c = map_files(argc, argv)) < 0) {
      fprintf(stderr, "Error mapping files\n");
      exit(0);
  }


  int i;

	/* destroy kseqs */
  for (i = 0; i < argc-4; i++) {

  	kseq_destroy(sequences[i]);
  }

	/* Free up */
	free(seq_l);
  kseq_destroy(map_sequence);

  free_hash_table(ht);
  free(minimizer_target_set.set);

}



int read_targets(int argc, char **argv) {

  if (argc < 4) {
    fprintf(stderr, "Need at least four arguments; <Targets indexed> <Queries> k w \n");
    exit(0);
  }

  /*
  * Currently, First n-1 files are the ones to index, last one is a file with query(ies)
  * Ex. "./minimap <index_file1> <index_file2> <query_file>
  */

  FILE *index_fp[argc-4];

  int i;
  int j;

	sequences = malloc((argc-4) * sizeof(kseq_t*));
	seq_l = malloc((argc-4) * sizeof(int));

	if (sequences == NULL) {

			fprintf(stderr, "Error allocating sequences\n");
			exit(0);
	}

  /* Read index files */

  for (i = 0; i < argc-4; i++) {

    index_fp[i] = fopen(argv[i+1], "r");

    if (index_fp[i] == 0) {
      perror("fopen");
      exit(1);
    }

    sequences[i] = kseq_init(fileno(index_fp[i]));

    while ((j = kseq_read(sequences[i])) >= 0) {
			seq_l[i] = sequences[i]->seq.l;
      #ifdef VERBOSE
	      printf("name: %s\n", sequences[i]->name.s);
				printf("len: %d\n", sequences[i]->seq.l);
      #endif VERBOSE
    }

    /* Close */
    fclose(index_fp[i]);
  }

  return 1;

}

int map_files(int argc, char ** argv) {

  /* Read map file */
  FILE *map_fp;
  map_fp = fopen(argv[argc-3], "r");
  int j;

  if (map_fp == 0) {
    perror("fopen");
    exit(1);
  }

  map_sequence = kseq_init(fileno(map_fp));

  while ((j = kseq_read(map_sequence)) >= 0) {
    #ifdef VERBOSE

    #endif VERBOSE

		/* Overlap -> map sequences */
    if (map_sequence->seq.l >= 2000 && map_sequence->seq.l <= 80000) {
		     map(map_sequence);
		     free(minimizer_query_set.set);
    }


  }

  /* Close */
  fclose(map_fp);

  return 1;
}

/* This is the lh3 function with support for my own code etc. */
int mm_sketch(char * s, minimizer_set_t * m_set, uint64_t w, uint64_t k, int sid) {

  int len = sid > 0 ? seq_l[sid] : strlen(s); //add kseq later.. TODO
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};

	int i, j, l, buf_pos, min_pos;
	minimizer_t *buf, min = { UINT64_MAX, UINT64_MAX };
  int set_counter = 0;

	assert(len > 0 && w > 0 && k > 0);
	buf = (minimizer_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
    int c = seq_nt4_table[(uint8_t)s[i]];
		minimizer_t info = { UINT64_MAX, UINT64_MAX };

    if (c < 4) { // not an ambiguous base

      int z;

      kmer[0] = (kmer[0] << 2 | c) & BIT_MASK;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer

      if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand

      if (++l >= k) {
          info.min_value = invertible_hash(kmer[z], BIT_MASK), info.idx = (uint64_t)sid<<32 | (uint32_t)i<<1 | z;
        }
		}

    else l = 0;

		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

    if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet

      for (j = buf_pos + 1; j < w; ++j) {

        if (min.min_value == buf[j].min_value && buf[j].idx != min.idx) {
            insert_minimizer(m_set, buf[j]);
            set_counter++;
        }
      }

			for (j = 0; j < buf_pos; ++j) {

        if (min.min_value == buf[j].min_value && buf[j].idx != min.idx) {

            insert_minimizer(m_set, buf[j]);
            set_counter++;
        }
      }
		}

		if (info.min_value <= min.min_value) { // a new minimum; then write the old min

      if (l >= w + k) {

        insert_minimizer(m_set, min);
        set_counter++;
      }
			min = info, min_pos = buf_pos;
		}
    else if (buf_pos == min_pos) { // old min has moved outside the window

    	if (l >= w + k - 1) {

          insert_minimizer(m_set, min);
          set_counter++;
      }

			for (j = buf_pos + 1, min.min_value = UINT64_MAX; j < w; ++j) { // the two loops are necessary when there are identical k-mers

           if (min.min_value >= buf[j].min_value) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
      }

			for (j = 0; j <= buf_pos; ++j) {

        if (min.min_value >= buf[j].min_value) min = buf[j], min_pos = j;
      }


			if (l >= w + k - 1) { // write identical k-mers

        for (j = buf_pos + 1; j < w; ++j) { // these two loops make sure the output is sorted

           if (min.min_value == buf[j].min_value && min.idx != buf[j].idx) {
              insert_minimizer(m_set, buf[j]);
              set_counter++;
           }
        }

				for (j = 0; j <= buf_pos; ++j) {
					if (min.min_value == buf[j].min_value && min.idx != buf[j].idx) {
            insert_minimizer(m_set, buf[j]);
            set_counter++;
          }
        }
			}
		}

    if (++buf_pos == w) {
      buf_pos = 0;
    }
	}
	if (min.min_value != UINT64_MAX) {
    insert_minimizer(m_set, min);
    set_counter++;
  }

  return set_counter;
}


/* Insert a minimizer */
inline void insert_minimizer(minimizer_set_t * m_set, minimizer_t v) {

	  if (m_set->ct == m_set->limit) {

				m_set->limit = (MINIMIZER_BLOCK * ++m_set->blk_ct);

	      if ((m_set->set = realloc(m_set->set, sizeof(minimizer_t) * MINIMIZER_BLOCK * (m_set->blk_ct))) == NULL) {
	        fprintf(stderr, "Error reallocating in minimizer_sketch\n");
	        exit(0);
	      }

	  }


    m_set->set[m_set->ct++] = v;

}

/* Initialises minimizer struct for sets used in indexing and mapping */
void init_minimizer_set(minimizer_set_t * m_set) {

    if ((m_set->set = malloc(sizeof(minimizer_t)*MINIMIZER_BLOCK)) == NULL) {

      fprintf(stderr, "Error allocating minimizer set in mm_index() A\n");
      exit(0);
    }

    m_set->ct = 0;
    m_set->blk_ct = 1;
		m_set->limit = MINIMIZER_BLOCK;

}

hash_table_t * mm_index() {

init_minimizer_set(&minimizer_target_set);

		/* Largest size of key */
		BIT_SIZE = opts.k << 1;
		BIT_MASK = ((uint64_t)(1 << BIT_SIZE)-1);
		const uint64_t KEY_SIZE = (1 << BIT_SIZE) - 1;

		int i;
		int cc = 0;

		for (i = 0; i < opts.n; i++) {
		  int c;

		  if ((c = mm_sketch(sequences[i]->seq.s, &minimizer_target_set, opts.w, opts.k, i)) < 0) {

		    fprintf(stderr, "Something went horribly wrong\n");
		    exit(0);
		  }

		  cc += c;

}


  if (cc < MINIMIZER_BLOCK*sizeof(minimizer_t)){
    if ((minimizer_target_set.set = realloc(minimizer_target_set.set, cc * sizeof(minimizer_t))) == NULL) {

        fprintf(stderr, "Error reallocating minimizer_target_set.set in index\n");
        exit(0);
    }

  }

  radix_sort_minimizer(minimizer_target_set.set, cc, 32/*BIT_SIZE*/);

  int64_t ct;

	/* Create hash table */
  ht = create_table(minimizer_target_set.ct);

  if (ht == NULL) {
    fprintf(stderr, "Error allocating hash table\n");
    exit(0);
  }

  uint64_t temp = minimizer_target_set.set[0].min_value;
	int temp_index = 0;
	int temp_count = 1;
  uint64_t v = NULL;


	/* Insert minimizers to hash table */
  for (ct = 1; ct < cc; ct++) {

      v = minimizer_target_set.set[ct].min_value;

			/* On new minimizers, insert index to previous (temp) */
      if (temp != v) {

        insert(ht, temp, temp_index, temp_count);
        temp = v;
				temp_index = ct;
				temp_count = 1;
      } else {

        temp_count++;
      }

  }
	/* Add last minimizer index */
	insert(ht, temp, temp_index, temp_count);

  return ht;

}


/* Prints the minimizer set */
void print_minimizer_table(minimizer_t * t, int len) {

  int count;

  for (count = 0; count < len; count++) {

    printf("%d: %llu-- seq %d idx %d strand %d\n", count, t[count].min_value, get_sequence(t[count].idx), get_index(t[count].idx), get_strand(t[count].idx));

  }
}

/* Prints overlap set of A_t */
void print_A_table(__uint128_t * A, int len) {

  int i;
	uint64_t idx;

  for (i = 0; i < len; i++) {

		if (get_A_strand(A[i]) == 0) {
			idx = get_A_diagonal(A[i]) - BIAS + get_A_index(A[i]);
		} else {
			idx = get_A_diagonal(A[i]) - BIAS - get_A_index(A[i]);
		}

    printf("%d -- seq: %u, strand: %u, diag: %d t_index: %u, q_indx: %u\n", i, get_A_sequence(A[i]), get_A_strand(A[i]), get_A_diagonal(A[i]) - BIAS, get_A_index(A[i]), idx);

  }
}

/* Prints an array of 128-bit values */
void print_128_array(__uint128_t * arr, int len) {

  int count;

  for (count = 0; count < len /* fix this later */; count++) {
    printf("%d ", count);
    print_uint128(arr[count]);
    printf("\n");

  }
}

/* Prints a 128-bit value with n-imal notation (2/10)*/
void print_uint128(__uint128_t n)
{

  if (n == 0) {
    return;
  }

  print_uint128(n/2);
  putchar(n%2+0x30);
}


/* overlap and map q to index */
void map(kseq_t * q) {

    __uint128_t * A;


  init_minimizer_set(&minimizer_query_set);

  int A_counter = 0;
  int A_blk_ct = 1;

  int cc;

  if ((cc = mm_sketch(q->seq.s, &minimizer_query_set, opts.w, opts.k, -1)) < 0) {
    fprintf(stderr, "Something went horribly wrong\n");
    exit(0);
  }

  if (cc == 0) {
      /* No minimizers */
			free(A);
      return;
  }

  /* Allocate to same size as query M array should be fine, as it cant grow any larger (???) */
  if ((A = malloc(sizeof(__uint128_t)*MINIMIZER_BLOCK*A_blk_ct)) == NULL) {

    fprintf(stderr, "Error allocating A in map()\n");
    exit(0);
  }

  int i;
  entry_t * entry;
  int max_occ = 8;

  /* Collect minimizer hits */
  for (i = 0; i < cc; i++) {

    int j;
    minimizer_t * m_q; // minimizer from query
    minimizer_t * m_i; // minimizer from index

    m_q = &(minimizer_query_set.set[i]);

    /* If theres a hit */
    if ((entry = get(ht, m_q->min_value)) != NULL) {

      uint64_t q_i;
      int32_t diag;

      /* Iterate over hits */
      for (j = 0; j < entry->count && j < max_occ; j++) {


        /* get index to the minimizer_target_set */
        q_i = entry->val.index + j;

        /* get minimizer based on index found in H[h] */
        m_i = &(minimizer_target_set.set[q_i]);

        if (A_counter == (MINIMIZER_BLOCK * A_blk_ct)) {

            if ((A = realloc(A, MINIMIZER_BLOCK * sizeof(__uint128_t) * ++A_blk_ct)) == NULL) {
                fprintf(stderr, "Error reallocating new A in map()\n");
                exit(0);
            }
        }

				/* Insert to overlap set */
        if (get_strand(m_q->idx) == get_strand(m_i->idx)) {

          diag = (get_index(m_q->idx) - get_index(m_i->idx)) + BIAS;

          A[A_counter++] = (((__uint128_t)get_sequence(m_i->idx)) << 96) | (/*(__uint128_t)get_strand(m_q->idx)*/((__uint128_t)0x0) << 64) | (((uint64_t)diag) << 32 | ((__uint128_t)get_index(m_i->idx))); /* (get_strand) or m_i, doesn't matter.. */

        } else {

          diag = (get_index(m_q->idx) + get_index(m_i->idx)) + BIAS;

          A[A_counter++] = (((__uint128_t)get_sequence(m_i->idx)) << 96) | (/*(__uint128_t)get_strand(m_q->idx)*/((__uint128_t)0x1) << 64) | (((uint64_t)diag) << 32 | ((__uint128_t)get_index(m_i->idx))); /* (get_strand) or m_i, doesn't matter.. */

        }

      }

    }

}

  if (A_counter == 0) {
			fprintf(stderr, "No matches with query %s\n", q->name.s);
      free(A);
      return;
  }

  if ((A = realloc(A, sizeof(__uint128_t)*A_counter)) == NULL) {
		printf("length: %d \n", A_counter);
    fprintf(stderr, "Error reallocating A\n");

    exit(0);
  }

  /* Sort A on sequence, strand, diagonal, target index. */
  radix_sort_128(A, A_counter, 128);

	map_overlaps(A, A_counter, q);

}

int map_overlaps (__uint128_t * A, int A_counter, kseq_t * q) {

	int b = 0;
	int e;

	int ct;

	int diff = 0;
	uint32_t * lis = malloc(1);
	uint32_t * p = malloc(1);

	/* Clusters within radius */
	for (e = 1; e < A_counter; e++) {
		// break; /* Uncomment for no map. ie. to compare to parallel.

		if ((e == A_counter - 1)
		|| (get_A_sequence(A[e + 1]) != get_A_sequence(A[e]))
		|| (get_A_strand(A[e + 1]) != get_A_strand(A[e]))
		|| (abs(get_A_diagonal(A[e + 1]) - get_A_diagonal(A[e]))) >= opts.r)
		{


			diff = e - b + 1; // +1 for 0-index
			// printf("b %d e %d diff %d\n", b, e, diff);
			/* Sort this interval on query index ... Qsort or radix?. Radix may still be faster
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

				/* sort on queries */
				qsort(&A[b], diff, sizeof(A_t), compare_query_index);

				/* start init LIS routine */

				if ((lis = realloc(lis, sizeof(uint32_t)*diff)) == NULL) {
						fprintf(stderr, "Error allocating lis in minimap.c\n");
						exit(0);
				}

				if ((p = realloc(p, sizeof(uint32_t)*diff)) == NULL) {
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

											print_chain(&A[b], lis, q, inner_b, inner_e);
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

							print_chain(&A[b], lis, q, inner_b, inner_e);
					}
				}

			}

			b = e + 1;


		}
	}

	free(lis);
	free(p);
	free(A);

}

/* Return the query index back from diagonal */
inline uint32_t convert_index(A_t a) {

  if (get_A_strand(a) == 0) {

    return get_A_diagonal(a) + get_A_index(a);

  } else {

    return get_A_diagonal(a) - get_A_index(a);
  }
}

/* Used with qsort in map() */
inline int compare_query_index(A_t * a, A_t * b) {

		return convert_index(*a) - convert_index(*b);
}

inline void print_chain(A_t * A, uint32_t * lis, kseq_t * q, int b, int e) {
	int sequence = get_A_sequence(A[b]);																																																																																				 // Need a fix for this...
	printf("%-25s %-15d %-9d %-9d %-3c %-25s %-9u %-9u \n", /*fix this*/q->name.s, q->seq.l, convert_index(A[lis[b]]) - BIAS - opts.k + 1, convert_index(A[lis[e]]) - BIAS + 1, "+-"[get_A_strand(A[b])], sequences[sequence]->name.s, get_A_index(A[lis[b]]) - opts.k + 1, get_A_index(A[lis[e]]) + 1);

}


/* For minimizer_t */

inline uint32_t get_index(uint64_t i) {

  return (uint32_t)((i & 0xFFFFFFFF) >> 1);
}

inline uint32_t get_strand(uint64_t i) {

  return (uint32_t)(i & 0x1);
}

inline uint32_t get_sequence(uint64_t i) {

  return (uint32_t)(i >> 32);
}

/* For A_t */

inline uint32_t get_A_sequence(__uint128_t i) {

  return (uint32_t)((i >> 96) & 0xFFFFFFFF);
}

inline uint32_t get_A_strand(__uint128_t i) {

  return (uint32_t)((i >> 64) & 0xFFFFFFFF);
}

inline uint32_t get_A_diagonal(__uint128_t i) {

  return (uint32_t)((i >> 32) & 0xFFFFFFFF);
}

inline uint32_t get_A_index(__uint128_t i) {

  return (uint32_t)(i & 0xFFFFFFFF);
}
