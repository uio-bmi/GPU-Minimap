#include <stdio.h>
#include <assert.h>
#include "minimap.h"

inline uint32_t less_than(uint32_t x, uint32_t y) {

    return x < y;
}

inline uint32_t greater_than(uint32_t x, uint32_t y) {

    return x > y;
}


/*
 * Finds the longest increasing subsequence of A, placing the subsequence's inde
 * -xes in an array s. strand determines whether the algorithm should sort by
 * < or >.
 *
 * slightly modified version from lh3/minimap/ksort.h, which in turn is translated from:
 * http://www.algorithmist.com/index.php/Longest_Increasing_Subsequence.cpp
 *
 */
int longest_increasing_subsequence(A_t * a, uint32_t * b, uint32_t * _p, int strand, uint32_t length) {

    uint32_t u, v, i;
    uint32_t * top = b;
    uint32_t * p = _p;

    assert(b != NULL);


    assert(length != 0);

    uint32_t (*func) (uint32_t x, uint32_t y) = strand ? greater_than : less_than;

    *top++ = 0;

    for (i = 1; i < length; i++) {

        if (func(get_A_index(a[*(top-1)]), get_A_index(a[i]))) {
            p[i] = *(top-1);
            *top++ = i;
            continue;
        }


        for (u = 0, v = top - b - 1; u < v;) {

            uint32_t c = (u + v) >> 1;

            if (func(get_A_index(a[b[c]]), get_A_index(a[i]))) {
                u = c + 1;
            }

            else {
                v = c;
            }
        }
        if (func(get_A_index(a[i]), get_A_index(a[b[u]]))) {
            if (u > 0) p[i] = b[u-1];
            b[u] = i;
        }

    }

    for (u = top - b, v = *(top-1); u--; v = p[v]) b[u] = v;

    return top - b;
}
