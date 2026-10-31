#ifndef ZEROIZE_H
#define ZEROIZE_H

#include <stddef.h>

/*
 * Best-effort secure zeroization that the compiler cannot remove as a
 * dead store: a volatile-qualified store loop is outside the reach of
 * dead-store elimination and has no libc dependency, so it stays
 * C99-portable across the reference and optimized trees.
 *
 * Used to erase secret key material (entropy seed, kappa, s = (s0,s1,e),
 * masking nonce y, c*s products, and the XOF states they feed) when it
 * goes out of scope.  See the commit message for the coverage map.
 */
static inline void darts_zeroize(void *p, size_t len)
{
    volatile unsigned char *vp = (volatile unsigned char *)p;
    size_t i;
    for (i = 0; i < len; i++) {
        vp[i] = 0;
    }
}

#endif /* ZEROIZE_H */
