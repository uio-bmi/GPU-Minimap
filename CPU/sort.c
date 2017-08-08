
#include <stdio.h>
#include <assert.h>
#include "minimap.h"


/*
 * TODO: Find an elegant way of generalizing the sort for both, as they have a lot of overlapping code.
 * val's really only need 8 bit types if we're gonna stick to 8-bit key sizes.
 *
 */

/* Radix sort implementation */

/* Sorts array of minimizers */
void radix_sort_minimizer(minimizer_t * src, uint64_t size, uint64_t max) {

    /* For switching between passes */
    minimizer_t * arr = malloc(size*sizeof(minimizer_t));

    minimizer_t * ptrs[2] = {src, arr}; //src = 0, arr = 1;

    /* This looks somewhat tedious, but it's mostly for keeping track of parameters atm.. */
    int key_size = 8; //max >> 2; //tweak later...
    int passes = max / key_size; // We will be sorting up to 64-bits (k <= 32) for now..Check back if upgrade to 128-bit
    int b_size = 1 << key_size; //256
    int b_mask = b_size - 1;

    /* Initialize buckets */
    int buckets[b_size];
    memset(&buckets, 0x0, sizeof(int)*b_size);

    /* For 64 bits, this makes it 8 passes. */


    int i;
    int j;


    // First pass always begins with 0
    int ptr_idx = 0;

    /* On 'jth pass (1st being 0...)*/
    for (i = 0; i < passes; i++) {

        /* Place values into buckets */
        for (j = 0; j < size; j++) {

            uint64_t val = (ptrs[ptr_idx][j].min_value >> i*key_size) & b_mask;
            buckets[val]++;

        }

        /* Align buckets */
        int c;

        for (c = 1; c < b_size; c++) {

            buckets[c] += buckets[c-1];
        }


        for (j = size-1; j >= 0; j--) {

            minimizer_t * temp = &ptrs[ptr_idx][j];
            uint64_t val = (temp->min_value >> i*key_size) & b_mask;

            ptrs[ptr_idx ^ 0x1][--buckets[val]] = *temp;
        }

        /* Flip index */
        ptr_idx ^= 0x1;

        memset(&buckets, 0x0, sizeof(int)*b_size);


    }

    free(arr);
}

/*
 * Sorts array of (128-bit) .. Will merge with specialised "minimizer_t" version
 * at some point...
 */
void radix_sort_128(__uint128_t * src, uint64_t size, uint64_t max) {

    /* For switching between passes */
    __uint128_t * arr = malloc(sizeof(__uint128_t)*size);

    __uint128_t * ptrs[2] = {src, arr}; //src = 0, arr = 1;

    /* we assume 128-bit max here... */
    int key_size = 8;//0x10; //tweak later...
    int passes = max / key_size; // We will be sorting up to 128-bits ..
    int b_size = 1 << key_size; //256
    int b_mask = b_size - 1;

    /* Initialize buckets */
    int * buckets;
    if ((buckets = calloc(b_size, sizeof(int))) == NULL) {
        fprintf(stderr, "error allocating buckets in radix_sort_128\n");
        exit(0);
    }

    /* For 64 bits, this makes it 8 passes. */

    int i;
    int j;

    // First pass always begins with 0
    int ptr_idx = 0;

    uint32_t b_val;
    __uint128_t * temp;
    __uint128_t val;


    /* On 'jth pass (1st being 0...)*/
    for (i = 0; i < passes; i++) {


        /* Place values into buckets */
        for (j = 0; j < size; j++) {

            b_val = (ptrs[ptr_idx][j] >> i*key_size) & b_mask;

            buckets[b_val]++;

        }

        /* Align buckets */
        int c;

        for (c = 1; c < b_size; c++) {

            buckets[c] += buckets[c-1];
        }

        for (j = size-1; j >= 0; j--) {

            temp = &ptrs[ptr_idx][j];

            val = (*temp >> i*key_size) & b_mask;

            ptrs[ptr_idx ^ 0x1][--buckets[val]] = *temp;
        }

        /* Flip index */
        ptr_idx ^= 0x1;

        memset(buckets, 0x0, sizeof(int)*b_size);


    }
    free(buckets);
    free(arr);
}
