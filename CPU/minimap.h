#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "kseq.h"

KSEQ_INIT(FILE*, read)


typedef struct {

  uint64_t min_value; /* minimizer, h */
  uint64_t idx;

} minimizer_t;

typedef struct {

    uint32_t * s;
    int size;
    int top;
    int count;

} vector_t;

typedef struct {

    minimizer_t * set;
    uint64_t ct;
    uint64_t blk_ct;
    uint64_t limit;

} minimizer_set_t;

typedef struct {

    int k; // K-mer
    int w; // window size
    int r; // (r) How large difference in diagonals allowed (500 default)
    float m; // minimizeroverlap % (default 0.50 (50%)); What does this mean ?: "Minimap merges two windows if -m=50% of minimizer matches overlap."
    int c; // minimum matches (chain length)
    int L; //matching length
    int f; // n/a
    int g; //Break the chain whenever there is a gap longer than -g=10000.
    int n; // Number of inputs (targets)

} options_t;

typedef struct {

  /* Value */
  uint64_t min_value;
  uint64_t index;

} value_t;

typedef struct entry_t entry_t;

struct entry_t {

  /* Key */
  uint64_t key;

  /* Value */
  value_t val;

  /* Count */
  uint64_t count;


};

/* The hash table */
typedef struct {

  int size;
  entry_t **table;

} hash_table_t;


/* A in map. contains t=target, r=strand, c=diagonal, i=index */
typedef __uint128_t A_t;

/* I/O */
int read_targets(int argc, char ** argv);
int map_files(int argc, char ** argv);

/* minimap.c stuff */
int minimizer_sketch(char * s, minimizer_set_t * set, uint64_t w, uint64_t k, int sid);
int mm_sketch(char * s, minimizer_set_t * m_set, uint64_t w, uint64_t k, int sid);
int compare_minimizer(minimizer_t a, minimizer_t b);
inline void insert_minimizer(minimizer_set_t * m_set, minimizer_t v);
void init_minimizer_set(minimizer_set_t * m_set);

hash_table_t * mm_index();
inline uint64_t invertible_hash(uint64_t x, uint64_t m);
uint64_t invertible_hashi(uint64_t x, uint64_t p);

void print_minimizer_table(minimizer_t * t, int len);
void print_128_array(__uint128_t * arr, int len);
void print_A_table(__uint128_t * A, int len);

void map(kseq_t * q);
int map_overlaps (__uint128_t * A, int A_counter, kseq_t * q);

inline uint32_t convert_index(A_t a);
int compare_query_index(A_t * a, A_t * b);
void print_chain(A_t * A, uint32_t * lis, kseq_t * q, int b, int e);

/* For minimizer_t */
inline uint32_t get_index(uint64_t i);
inline uint32_t get_strand(uint64_t i);
inline uint32_t get_sequence(uint64_t i);

/* For 128 bit (A_t) */
inline uint32_t get_A_sequence(__uint128_t i);
inline uint32_t get_A_strand(__uint128_t i);
inline uint32_t get_A_diagonal(__uint128_t i);
inline uint32_t get_A_index(__uint128_t i);

/* misc. */
void print_uint128(__uint128_t n);

/* hash.c stuff */
hash_table_t * create_table(int size);
void insert(hash_table_t * ht, uint64_t key, uint64_t idx, int count);
inline entry_t * get(hash_table_t * ht, uint64_t key);
void print_with_key(uint64_t key, hash_table_t * ht, minimizer_t * m);
uint64_t calc_needed_bits(uint64_t integer);
uint64_t get_lookup_key(hash_table_t * ht, uint64_t key);
void free_hash_table(hash_table_t * ht);

/* sort.c stuff */
void radix_sort_minimizer(minimizer_t * src, uint64_t size, uint64_t max);
void radix_sort_128(__uint128_t * src, uint64_t size, uint64_t max);

/* LIS.c */
uint32_t less_than(uint32_t x, uint32_t y);
uint32_t greater_than(uint32_t x, uint32_t y);

void init_vector(vector_t * vector, uint32_t size);
uint32_t back(vector_t * vector);
void push(vector_t * vector, uint32_t val);
void vector_insert(vector_t * vec, uint32_t i, uint32_t val);
int longest_increasing_subsequence(A_t * a, uint32_t * b, uint32_t * p, int strand, uint32_t length);
