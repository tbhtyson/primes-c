#include "define.h"
#include <stdlib.h>
#include <signal.h>
// Miller-Rabin primality test
// more rounds = more accurate; 20 rounds is very reliable
int is_prime_miller_rabin(const uint256_t *n, int rounds) {
    uint256_t *two  = calloc(1, sizeof(uint256_t));
    uint256_t *one  = calloc(1, sizeof(uint256_t));
    uint256_t *d    = calloc(1, sizeof(uint256_t));
    uint256_t *nm1  = calloc(1, sizeof(uint256_t));
    uint256_t *a    = calloc(1, sizeof(uint256_t));
    uint256_t *x    = calloc(1, sizeof(uint256_t));
    uint256_t *q    = calloc(1, sizeof(uint256_t));
    uint256_t *r    = calloc(1, sizeof(uint256_t));

    one->limbs[0] = 1;
    two->limbs[0] = 2;

    int result = 1;

    if (cmp256(n, two) < 0)  { result = 0; goto done; }
    if (cmp256(n, two) == 0) { result = 1; goto done; }
    if (n->limbs[0] % 2 == 0){ result = 0; goto done; }

    // write n-1 as 2^s * d
    sub256(n, one, nm1);
    *d = *nm1;
    uint64_t s = 0;
    while (d->limbs[0] % 2 == 0) {
        shift_right256(d, 1, d);
        s++;
    }

    // test with small fixed witnesses (good enough for 256-bit)
    uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                             31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int num_witnesses = rounds < 20 ? rounds : 20;

    for (int w = 0; w < num_witnesses; w++) {
        a->limbs[0] = witnesses[w];
        if (cmp256(a, n) >= 0) continue;

        powmod256(a, d, n, x);

        if (cmp256(x, one) == 0 || cmp256(x, nm1) == 0) continue;

        int composite = 1;
        for (uint64_t j = 0; j < s - 1; j++) {
            mulmod256(x, x, n, r);
            *x = *r;
            if (cmp256(x, nm1) == 0) { composite = 0; break; }
        }

        if (composite) { result = 0; goto done; }
    }

done:
    free(two); free(one); free(d); free(nm1);
    free(a); free(x); free(q); free(r);
    return result;
}

volatile sig_atomic_t keep_running = 1;

void handle_sigint(int sig) {
    printf("\nCaught signal %d (Ctrl+C). Cleaning up...\n", sig);
    keep_running = 0;
}

int main() {
  uint256_t *n = calloc(1, sizeof(uint256_t));
  uint256_t *o = calloc(1, sizeof(uint256_t));
  uint256_t *p = calloc(1, sizeof(uint256_t));
  (*n).limbs[0] = 3;
  (*o).limbs[0] = 2;
  signal(SIGINT, handle_sigint);
  while(keep_running) {
    if(is_prime_miller_rabin(n, 20)) {
      bigprint(n);
    }
    add256(n, o, p);
    *n = *p;
  }
  free(n);
  free(o);
  free(p);
}
