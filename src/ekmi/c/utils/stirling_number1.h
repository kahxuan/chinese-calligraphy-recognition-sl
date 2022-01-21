// https://rosettacode.org/wiki/Stirling_numbers_of_the_first_kind#C

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
 
typedef struct stirling_cache_tag1 {
    int max;
    int* values;
} stirling_cache1;
 
int stirling_number1(stirling_cache1* sc, int n, int k) {
    if (k == 0)
        return n == 0 ? 1 : 0;
    if (n > sc->max || k > n)
        return 0;
    return sc->values[n*(n-1)/2 + k - 1];
}
 
bool stirling_cache_create1(stirling_cache1* sc, int max) {
    int* values = calloc(max * (max + 1)/2, sizeof(int));
    if (values == NULL)
        return false;
    sc->max = max;
    sc->values = values;
    for (int n = 1; n <= max; ++n) {
        for (int k = 1; k <= n; ++k) {
            int s1 = stirling_number1(sc, n - 1, k - 1);
            int s2 = stirling_number1(sc, n - 1, k);
            values[n*(n-1)/2 + k - 1] = s1 + s2 * (n - 1);
        }
    }
    return true;
}
 
void stirling_cache_destroy1(stirling_cache1* sc) {
    free(sc->values);
    sc->values = NULL;
}


