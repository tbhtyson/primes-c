#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    uint64_t limbs[4];   /* limbs[0] = least significant */
} uint256_t;

int is_zero256(const uint256_t *a) {
    for (int i = 0; i < 4; i++)
        if (a->limbs[i]) return 0;
    return 1;
}
 
/* returns 1 if a > b, -1 if a < b, 0 if equal */
int cmp256(const uint256_t *a, const uint256_t *b) {
    for (int i = 3; i >= 0; i--) {
        if (a->limbs[i] > b->limbs[i]) return  1;
        if (a->limbs[i] < b->limbs[i]) return -1;
    }
    return 0;
}
 
 
/* ─────────────────────────────────────────────
   Shift (internal tmp avoids aliasing bugs)
   ───────────────────────────────────────────── */
 
void shift_left256(const uint256_t *a, int bits, uint256_t *result) {
    uint256_t tmp = {{0}};
    int limb_shift = bits / 64;
    int bit_shift  = bits % 64;
 
    for (int i = 3; i >= 0; i--) {
        if (i >= limb_shift)
            tmp.limbs[i] = a->limbs[i - limb_shift] << bit_shift;
        if (bit_shift && i > limb_shift)
            tmp.limbs[i] |= a->limbs[i - limb_shift - 1] >> (64 - bit_shift);
    }
    *result = tmp;
}
 
void shift_right256(const uint256_t *a, int bits, uint256_t *result) {
    uint256_t tmp = {{0}};
    int limb_shift = bits / 64;
    int bit_shift  = bits % 64;
 
    for (int i = 0; i < 4; i++) {
        if (i + limb_shift < 4)
            tmp.limbs[i] = a->limbs[i + limb_shift] >> bit_shift;
        if (bit_shift && i + limb_shift + 1 < 4)
            tmp.limbs[i] |= a->limbs[i + limb_shift + 1] << (64 - bit_shift);
    }
    *result = tmp;
}
 
 
/* ─────────────────────────────────────────────
   Addition / Subtraction
   ───────────────────────────────────────────── */
 
void add256(const uint256_t *a, const uint256_t *b, uint256_t *result) {
    uint64_t carry = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t ai = a->limbs[i];
        uint64_t bi = b->limbs[i];
        result->limbs[i] = ai + bi + carry;
        carry = (result->limbs[i] < ai) || (carry && result->limbs[i] == ai);
    }
}
 
void sub256(const uint256_t *a, const uint256_t *b, uint256_t *result) {
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t ai = a->limbs[i];
        uint64_t bi = b->limbs[i];
        result->limbs[i] = ai - bi - borrow;
        borrow = (ai < bi + borrow) ? 1 : 0;
    }
}
 
 
/* ─────────────────────────────────────────────
   Multiplication
   ───────────────────────────────────────────── */
 
/* multiply two 64-bit values without __uint128_t
   returns low 64 bits, stores high 64 bits in *hi */
static uint64_t mul64(uint64_t a, uint64_t b, uint64_t *hi) {
    uint32_t a_lo = (uint32_t)a,        a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)b,        b_hi = (uint32_t)(b >> 32);
 
    uint64_t ll = (uint64_t)a_lo * b_lo;
    uint64_t lh = (uint64_t)a_lo * b_hi;
    uint64_t hl = (uint64_t)a_hi * b_lo;
    uint64_t hh = (uint64_t)a_hi * b_hi;
 
    uint64_t mid   = lh + hl;
    uint64_t carry = (mid < lh) ? 1 : 0;
    uint64_t lo    = ll + (mid << 32);
    carry += (lo < ll) ? 1 : 0;
 
    *hi = hh + (mid >> 32) + (carry << 32);
    return lo;
}
 
void mul256(const uint256_t *a, const uint256_t *b, uint256_t *result) {
    for (int i = 0; i < 4; i++) result->limbs[i] = 0;
 
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4 - i; j++) {
            uint64_t hi, lo = mul64(a->limbs[i], b->limbs[j], &hi);
            uint64_t prev = result->limbs[i+j];
            result->limbs[i+j] += lo + carry;
            carry = hi + (result->limbs[i+j] < prev ? 1 : 0);
        }
    }
}
 
 
/* ─────────────────────────────────────────────
   Division
   ───────────────────────────────────────────── */
 
void div256(const uint256_t *a, const uint256_t *b, uint256_t *q, uint256_t *r) {
    for (int i = 0; i < 4; i++) q->limbs[i] = r->limbs[i] = 0;
 
    for (int i = 255; i >= 0; i--) {
        uint256_t shifted;
        shift_left256(r, 1, &shifted);
        *r = shifted;
 
        int limb = i / 64, bit = i % 64;
        if (a->limbs[limb] & ((uint64_t)1 << bit))
            r->limbs[0] |= 1;
 
        if (cmp256(r, b) >= 0) {
            sub256(r, b, &shifted);
            *r = shifted;
            q->limbs[limb] |= (uint64_t)1 << bit;
        }
    }
}
 
 
/* ─────────────────────────────────────────────
   Square root (Newton's method)
   ───────────────────────────────────────────── */
 
void sqrt256(const uint256_t *n, uint256_t *result) {
    if (is_zero256(n)) { *result = *n; return; }
 
    uint256_t x = {{0}}, x2 = {{0}}, q = {{0}}, r = {{0}};
    shift_right256(n, 128, &x);
    if (is_zero256(&x)) x.limbs[0] = 1;
 
    while (1) {
        div256(n, &x, &q, &r);
        add256(&x, &q, &x2);
        shift_right256(&x2, 1, &x2);
        if (cmp256(&x2, &x) >= 0) break;
        x = x2;
    }
    *result = x;
}
 
 
/* ─────────────────────────────────────────────
   Modular arithmetic (no div256 in hot path)
   ───────────────────────────────────────────── */
 
void mulmod256(const uint256_t *a, const uint256_t *b, const uint256_t *m, uint256_t *result, uint256_t *ta, uint256_t *tb, uint256_t *tmp) {
 
    /* reduce a mod m by subtraction */
    *ta = *a;
    while (cmp256(ta, m) >= 0)
        sub256(ta, m, ta);
 
    *tb = *b;
    for (int i = 0; i < 4; i++) result->limbs[i] = 0;
 
    while (!is_zero256(tb)) {
        if (tb->limbs[0] & 1) {
            add256(result, ta, result);
            if (cmp256(result, m) >= 0)
                sub256(result, m, result);
        }
        shift_left256(ta, 1, tmp);
        *ta = *tmp;
        if (cmp256(ta, m) >= 0)
            sub256(ta, m, ta);
        shift_right256(tb, 1, tmp);
        *tb = *tmp;
    }
 
}
 
void powmod256(const uint256_t *base, const uint256_t *exp, const uint256_t *mod, uint256_t *result, uint256_t *b, uint256_t *e, uint256_t *tmp, uint256_t *ta, uint256_t *tb,     uint256_t *tmp2) {
    *e = *exp;
    result->limbs[0] = 1;
    result->limbs[1] = result->limbs[2] = result->limbs[3] = 0;
 
    /* reduce base mod m by subtraction */
    *b = *base;
    while (cmp256(b, mod) >= 0)
        sub256(b, mod, b);
 
    while (!is_zero256(e)) {
        if (e->limbs[0] & 1) {
            mulmod256(result, b, mod, tmp, ta, tb, tmp2);
            *result = *tmp;
        }
        mulmod256(b, b, mod, tmp, ta, tb, tmp2);
        *b = *tmp;
        shift_right256(e, 1, tmp);   /* e >>= 1, no div256 needed */
        *e = *tmp;
    }
}

void print256_decimal(const uint256_t *n) {
    char buf[79];
    int  pos = 78;
    buf[78] = '\0';
 
    if (is_zero256(n)) { printf("0\n"); return; }
 
    uint256_t *cur = (uint256_t *)calloc(1, sizeof(uint256_t));
    uint256_t *ten = (uint256_t *)calloc(1, sizeof(uint256_t));
    uint256_t *q   = (uint256_t *)calloc(1, sizeof(uint256_t));
    uint256_t *r   = (uint256_t *)calloc(1, sizeof(uint256_t));
 
    *cur = *n;
    ten->limbs[0] = 10;
 
    while (!is_zero256(cur)) {
        div256(cur, ten, q, r);
        buf[--pos] = '0' + (char)r->limbs[0];
        *cur = *q;
    }
 
    printf("%s\n", &buf[pos]);
    free(cur); free(ten); free(q); free(r);
}


#endif
/* RULES: 
INITIALIZE ALL UINT256_T VALUES TO 0 BEFORE USE, OTHERWISE SLOP STARTS!
 */
