#include "define.h"
#include <stdlib.h>
#include <signal.h>

uint64_t mod256_small(const uint256_t *n, uint64_t m) {
    uint64_t r = 0;
    for (int i = 3; i >= 0; i--) {
        // process each 64-bit limb using 128-bit intermediate
        // split into two 32-bit halves to avoid overflow
        uint32_t hi = (uint32_t)(n->limbs[i] >> 32);
        uint32_t lo = (uint32_t)(n->limbs[i]);
        r = ((r * (0x100000000ULL % m) + hi) % m
                * (0x100000000ULL % m) + lo) % m;
    }
    return r;
}

void mod256_batch(const uint256_t *n, const uint64_t *primes, 
                  uint64_t *remainders, int count) {
    for (int p = 0; p < count; p++) remainders[p] = 0;

    for (int i = 3; i >= 0; i--) {
        uint32_t hi = (uint32_t)(n->limbs[i] >> 32);
        uint32_t lo = (uint32_t)(n->limbs[i]);
        for (int p = 0; p < count; p++) {
            uint64_t m = primes[p];
            remainders[p] = ((remainders[p] * (0x100000000ULL % m) + hi) % m
                            * (0x100000000ULL % m) + lo) % m;
        }
    }
}

static const uint64_t small_primes[] = {
    5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,181,191,193,197,199,211,223,227,
    229,233,239,241,251,257,263,269,271,277,281,283,293,307,
    311,313,317,331,337,347,349,353,359,367,373,379,383,389,
    397,401,409,419,421,431,433,439,443,449,457,461,463,467,
    479,487,491,499,503,509,521,523,541
};

int has_small_factor(const uint256_t *n) {
    for (int i = 0; i < 20; i++)
        if (mod256_small(n, small_primes[i]) == 0 
            && n->limbs[0] != small_primes[i])
            return 1;
    return 0;
}


int has_small_factor_batch(const uint256_t *n) {
    uint64_t remainders[100];
    mod256_batch(n, small_primes, remainders, 100);
    for (int i = 0; i < 100; i++)
        if (remainders[i] == 0 && n->limbs[0] != small_primes[i])
            return 1;
    return 0;
}

// Miller-Rabin primality test
// more rounds = more accurate; 20 rounds is very reliable
int is_prime_miller_rabin(const uint256_t *n, int rounds, uint256_t *d, uint256_t *nm1, uint256_t *a, uint256_t *x, uint256_t *r, uint256_t *b, uint256_t *e, uint256_t *tmp, uint256_t *ta, uint256_t *tb, uint256_t *tmp2, uint256_t *one, uint256_t *two) {
  for(int lx = 0; lx < 4; lx++) {
    (*a).limbs[lx]    = 0;
    /* tmp vars for the functions so calloc isn't slowing everything down by fucking 400% every time i run powmod and mulmod */
    (*b).limbs[lx]    = 0;
    (*e).limbs[lx]    = 0;
    (*tmp).limbs[lx]  = 0;
    (*ta).limbs[lx]   = 0;
    (*tb).limbs[lx]   = 0;
    (*tmp2).limbs[lx] = 0;
  }
    int result = 1;

    if (cmp256(n, two) < 0)  { result = 0; goto done; }
    if (cmp256(n, two) == 0) { result = 1; goto done; }
    if (n->limbs[0] % 2 == 0){ result = 0; goto done; }

    // write n-1 as 2^s * d
    sub256(n, one, nm1);
    *d = *nm1;
    uint64_t s = 0;
    while (d->limbs[0] % 2 == 0) {
        shift_right256(d, 1, tmp);
        *d = *tmp;
        s++;
    }

    // test with small fixed witnesses (good enough for 256-bit)
    static const uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                             31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int num_witnesses = rounds < 20 ? rounds : 20;

    for (int w = 0; w < num_witnesses; w++) {
        a->limbs[0] = witnesses[w];
        a->limbs[1] = a->limbs[2] = a->limbs[3] = 0;
        if (cmp256(a, n) >= 0) continue;

        powmod256(a, d, n, x, b, e, tmp, ta, tb, tmp2);

        if (cmp256(x, one) == 0 || cmp256(x, nm1) == 0) continue;

        int composite = 1;
        for (uint64_t j = 0; j < s - 1; j++) {
            mulmod256(x, x, n, r, ta, tb, tmp);
            *x = *r;
            if (cmp256(x, nm1) == 0) { composite = 0; break; }
        }

        if (composite) { result = 0; goto done; }
    }

done:
    return result;
}

int needed_rounds(const uint256_t *n) {
    // if upper 3 limbs are zero, n fits in 64 bits
    if (!n->limbs[3] && !n->limbs[2] && !n->limbs[1])
        return 3;   // 3 witnesses sufficient for n < 3,215,031,751
    if (!n->limbs[3] && !n->limbs[2])
        return 6;   // 6 for 128-bit range
    return 12;      // 12 for full 256-bit range
}

volatile sig_atomic_t keep_running = 1;

void handle_sigint(int sig) {
    printf("\nCaught signal %d (Ctrl+C). Cleaning up...\n", sig);
    keep_running = 0;
}

int main() {
  uint256_t *n = calloc(1, sizeof(uint256_t));
  uint256_t *p = calloc(1, sizeof(uint256_t));
  uint256_t *step = calloc(1, sizeof(uint256_t));
  uint256_t *four = calloc(1, sizeof(uint256_t));
    uint256_t *two  = calloc(1, sizeof(uint256_t));
    uint256_t *one  = calloc(1, sizeof(uint256_t));
    uint256_t *d    = calloc(1, sizeof(uint256_t));
    uint256_t *nm1  = calloc(1, sizeof(uint256_t));
    uint256_t *a    = calloc(1, sizeof(uint256_t));
    uint256_t *x    = calloc(1, sizeof(uint256_t));
    uint256_t *r    = calloc(1, sizeof(uint256_t));
    /* tmp vars for the functions so calloc isn't slowing everything down by fucking 400% every time i run powmod and mulmod */
    uint256_t *b    = calloc(1, sizeof(uint256_t));
    uint256_t *e    = calloc(1, sizeof(uint256_t));
    uint256_t *tmp  = calloc(1, sizeof(uint256_t));
    uint256_t *ta   = calloc(1, sizeof(uint256_t));
    uint256_t *tb   = calloc(1, sizeof(uint256_t));
    uint256_t *tmp2 = calloc(1, sizeof(uint256_t));

    one->limbs[0] = 1;
    two->limbs[0] = 2;
    four->limbs[0] = 4;

  printf("2\n3\n5\n");
  static const uint64_t wheel30[] = {4,2,4,2,4,6,2,6};  // gaps between 30k+1 candidates
  int wheel_pos = 0;
  
  uint64_t remainders[100];
  mod256_batch(n, small_primes, remainders, 100);

  // start at 7 (first value after handling 2,3,5 manually)
  (*n).limbs[0] = 7;
  signal(SIGINT, handle_sigint);
  while(keep_running) {
    // check small factors
    int composite = 0;
    for (int i = 0; i < 100; i++) {
        if (remainders[i] == 0 && n->limbs[0] != small_primes[i]) {
            composite = 1;
            break;
        }
    }    
    if(!composite && is_prime_miller_rabin(n, needed_rounds(n), d, nm1, a, x, r, b, e, tmp, ta, tb, tmp2, one, two)) {
      print256_decimal(n);
    }

    uint64_t gap = wheel30[wheel_pos];
    step->limbs[0] = gap;
    step->limbs[1] = step->limbs[2] = step->limbs[3] = 0;
    add256(n, step, p);
    *n = *p;
    wheel_pos = (wheel_pos + 1) % 8;

    // update remainders by the same gap
    for (int i = 0; i < 100; i++)
        remainders[i] = (remainders[i] + gap) % small_primes[i];
  }
  free(p);
  free(d); free(nm1); free(a); free(x); free(r); free(b); free(e); free(tmp); free(ta); free(tb); free(tmp2); free(n); free(step); free(four); free(two); free(one);
  
}
