
#include <stdio.h>
#include <assert.h>
#include "minimap.h"

int MASK;
#define LOAD_FACTOR 0.75


/* Creates hash table */
hash_table_t * create_table(int size) {


  if (size < 2) return NULL;

  hash_table_t * new_table = NULL;
  size = size / LOAD_FACTOR;

  if ((new_table = malloc(sizeof(hash_table_t))) == NULL) {
    fprintf(stderr, "Error\n");
    return NULL;
  }

  if ((new_table->table = calloc(size, sizeof(entry_t*))) == NULL) {
    fprintf(stderr, "Error\n");
    return NULL;
  }

  new_table->size = size;

  MASK = size;


  return new_table;

}

/* Algorithm 2 */
uint64_t invertible_hash(uint64_t x, uint64_t m) {

  x = (~x + (x << 21)) & m;
  x = x ^ x >> 24;
  x = (x + (x << 3) + (x << 8)) & m;
  x = x ^ x >> 14;
  x = (x + (x << 2) + (x << 4)) & m;
  x = x ^ x >> 28;
  x = (x + (x << 31)) & m;
  return x;

}

// uint64_t invertible_hashi(uint64_t x, uint64_t p) {

//   uint64_t m = ((uint64_t)(1 << p)-1);

//   uint64_t tmp;

// 	// Invert key = key + (key << 31)
// 	tmp = (x - (x << 31));
// 	x = (x - (tmp << 31)) & m;

// 	// Invert x = x ^ (x >> 28)
// 	tmp = x ^ x >> 28;
// 	x = x ^ tmp >> 28;

// 	// Invert x *= 21
// 	x = (x * 14933078535860113213ull) & m;

// 	// Invert x = x ^ (x >> 14)
// 	tmp = x ^ x >> 14;
// 	tmp = x ^ tmp >> 14;
// 	tmp = x ^ tmp >> 14;
// 	x = x ^ tmp >> 14;

// 	// Invert x *= 265
// 	x = (x * 15244667743933553977ull) & m;

// 	// Invert x = x ^ (x >> 24)
// 	tmp = x ^ x >> 24;
// 	x = x ^ tmp >> 24;

// 	// Invert x = (~x) + (x << 21)
// 	tmp = ~x;
// 	tmp = ~(x - (tmp << 21));
// 	tmp = ~(x - (tmp << 21));
// 	x = ~(x - (tmp << 21)) & m;

// 	return x;

// }

inline uint64_t hash_function(uint64_t x, uint64_t m) {
    return x % m;
}

/* Key is hash value */
/* uses quadratic probing,  */
void insert(hash_table_t * ht, uint64_t key, uint64_t idx, int count) {

  int i;

  /* AKA hash */
  uint64_t lookup_key = hash_function(key, MASK);

  for (i = 0; i < ht->size; i++) {

    uint64_t temp_key = (lookup_key + (i*i)) % MASK; // Wrap around..

    /* Has this bin 'bin' (pun intended) used yet? */
    if (ht->table[temp_key] == NULL) {

      if ((ht->table[temp_key] = malloc(sizeof(entry_t))) == NULL) {
        fprintf(stderr, "Error allocating space..\n");
        exit(0);
      }

      ht->table[temp_key]->key = temp_key; // need for checking later..
      ht->table[temp_key]->val = (value_t) {key, idx}; 
      ht->table[temp_key]->count = count;


      return;

    }

    /* If collision */
    else {

      /* If same actual key, then just reset the index */
      if (ht->table[temp_key]->val.min_value == key) {

        // This shouldn't even happen :)))

        ht->table[temp_key]->val = (value_t) {key, idx};
        ht->table[temp_key]->count++;
        return;
      }


    }

  }


}

/* Get log 2 base of s */
uint64_t calc_needed_bits(uint64_t s) {

    uint64_t c = 1;

    while (s >>= 1) {
        c++;
    }

    return c;
}

inline entry_t * get(hash_table_t * ht, uint64_t key) {


  if (ht->table[hash_function(key,MASK)] == NULL) {
    // Empty table..
    return NULL;
  }

  int i;
  uint64_t temp_key;

  for (i = 0; i < ht->size; i++) {

    temp_key = (key + (i*i)) % MASK;

    entry_t * entry = ht->table[temp_key];

    /* need this ? */
    if (entry != NULL) {

        if (entry->val.min_value == key) {

          return entry;
        }

    } else {
        //assert(0);
        return NULL;
    }
  }

  return (entry_t*) NULL;
}

void print_with_key(uint64_t key, hash_table_t * ht, minimizer_t * m) {


  entry_t * e = get(ht, key);

  if (e == NULL) {

    printf("Entry with key: %llu doesn't exist in this hash table.\n", key);
    return;
  }

  uint64_t index = e->val.index;

  int ck = 1;
  int i = 0;

  /* TODO: Rework with count later */
  while (ck) {
    if (m[index+i].min_value == key) {

      printf("Key: %llu, Index: %llu Minimiser: %llu Value: %llu \n", key, index+i, m[index+i].min_value, m[index+i].idx);
      i++;

    } else {

      ck = 0;
    }

  }

}

void free_hash_table(hash_table_t * ht) {

    int i;

    for (i = 0; i < ht->size; i++) {

        free(ht->table[i]);
    }

    free(ht->table);
    free(ht);
}
