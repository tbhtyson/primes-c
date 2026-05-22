#include "define.h"
#include <stdlib.h>
#include <signal.h>
void mulmod256(const uint256_t *a, const uint256_t *b,
               const uint256_t *m, uint256_t *result) {
    uint256_t *ta  = calloc(1, sizeof(uint256_t));
    uint256_t *tb  = calloc(1, sizeof(uint256_t));
    uint256_t *tmp = calloc(1, sizeof(uint256_t));
    uint256_t *q   = calloc(1, sizeof(uint256_t));

    // reduce a mod m first so ta < m
    div256(a, m, q, ta);
    *tb = *b;
    for (int i = 0; i < 4; i++) result->limbs[i] = 0;

    while (!is_zero256(tb)) {
        if (tb->limbs[0] & 1) {             // if b is odd
            add256(result, ta, result);
            if (cmp256(result, m) >= 0)     // keep result < m
                sub256(result, m, result);
        }
        shift_left256(ta, 1, tmp);          // ta = ta * 2
        *ta = *tmp;
        if (cmp256(ta, m) >= 0)             // keep ta < m
            sub256(ta, m, ta);
        shift_right256(tb, 1, tmp);         // tb = tb / 2
        *tb = *tmp;
    }

    free(ta); free(tb); free(tmp); free(q);
}

// use shift instead of div to halve e - much faster
// temporarily add these to find where it hangs
void powmod256(const uint256_t *base, const uint256_t *exp,
               const uint256_t *mod, uint256_t *result) {
    uint256_t *b   = calloc(1, sizeof(uint256_t));
    uint256_t *e   = calloc(1, sizeof(uint256_t));
    uint256_t *tmp = calloc(1, sizeof(uint256_t));
    uint256_t *q   = calloc(1, sizeof(uint256_t));

    *e = *exp;
    result->limbs[0] = 1;
    result->limbs[1] = result->limbs[2] = result->limbs[3] = 0;

    div256(base, mod, q, b);

    int iter = 0;
    while (!is_zero256(e)) {
        if (e->limbs[0] & 1) {
            mulmod256(result, b, mod, tmp);
            *result = *tmp;
        }
        mulmod256(b, b, mod, tmp);
        *b = *tmp;
        shift_right256(e, 1, tmp);
        *e = *tmp;
    }

    free(b); free(e); free(tmp); free(q);
}// Miller-Rabin primality test
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
