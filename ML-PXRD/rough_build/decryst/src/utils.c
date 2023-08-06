#ifdef BENCHMARK
#define _POSIX_C_SOURCE 199309L  // `struct timespec'.
#endif

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <arpa/inet.h>
#ifdef BENCHMARK
#include <time.h>
#endif

#include "utils.h"
#ifdef BENCHMARK
#include "bench_utils.h"
#endif

uint32_t hftol (float x) {
	unsigned neg = signbit (x) != 0;
	if (isnan (x)) return neg << 31 | 0xFF << 23 | 0x1;
	if (isinf (x)) return neg << 31 | 0xFF << 23;
	signed expo;
	float frac = frexpf (neg ? -x : x, &expo);
	if (frac == 0.0) return neg << 31;
	expo += 0x7E;
	return neg << 31 | (expo > 0 ?
		expo << 23 | ((uint32_t) (frac * 0x01000000) & 0x007FFFFF) :
		(uint32_t) (frac * (0x01000000 >> (1 - expo)))
	);
}

float hltof (uint32_t x) {
	signed sig = x >> 31 ? -1 : 1, expo = x >> 23 & 0xFF;
	unsigned frac = x & 0x007FFFFF;
	if (expo == 0xFF) return copysign (frac ? NAN : INFINITY, sig);
	return sig * (expo ?
		ldexpf (1.1920928955078125e-07 * (frac | 0x00800000), expo - 0x7F) :
		ldexpf (1.1920928955078125e-07 * frac, -0x7E)
	);
}

uint32_t hftonl (float x) {
	return htonl (hftol (x));
}

float nltohf (uint32_t x) {
	return hltof (ntohl (x));
}

void hftonl2 (void *buf, unsigned n) {
	float *src = buf;
	uint32_t *dest = buf;
	for (unsigned i = 0; i < n; ++i) dest[i] = hftonl (src[i]);
}

void nltohf2 (void *buf, unsigned n) {
	float *dest = buf;
	uint32_t *src = buf;
	for (unsigned i = 0; i < n; ++i) dest[i] = nltohf (src[i]);
}

int arefinite (float const buf[], unsigned n) {
	for (unsigned i = 0; i < n; ++i) if (!isfinite (buf[i])) return 0;
	return 1;
}

void cursor_mk (struct cursor *curs, FILE *file) {
	curs->file = file;
	curs->idx = -1;
}

void cursor_err (struct cursor *curs, char const desc[]) {
	fprintf (
		stderr, "%s near %u,%u\n",
		desc, curs->idx, (unsigned) (curs->ptr - curs->buf)
	);
}

int cursor_step (struct cursor *curs) {
	while (1) {
		curs->idx++;
		curs->ptr = curs->buf;
		if (!fgets (curs->buf, LEN_MAX, curs->file)) return 0;
		curs->len = strlen (curs->buf);
		if (curs->len && curs->buf[curs->len - 1] == '\n') {
			curs->buf[--curs->len] = '\0';
		} else cursor_err (curs, "W: oversized or incomplete line");
		if (!curs->len || curs->buf[0] != '#') return 1;
	}
}

int scan_space (char const **ptr) {
	char const *cur = *ptr;
	int ret = 0;
	while (1) {
		char c = *cur;
		if (c == ' ' || c == '\t') {
			++cur;
			ret = 1;
		} else break;
	}
	*ptr = cur;
	return ret;
}

int scan_unsigned (char const **ptr, unsigned *x) {
	unsigned long y;
	char *cur;
	errno = 0;
	*x = y = strtoul (*ptr, &cur, 10);
	int ret = !(cur == *ptr || errno) && y <= UINT_MAX;
	*ptr = cur;
	return ret;
}

int scan_signed (char const **ptr, signed *x) {
	signed long y;
	char *cur;
	errno = 0;
	*x = y = strtol (*ptr, &cur, 10);
	int ret = !(cur == *ptr || errno) && labs (y) <= INT_MAX;
	*ptr = cur;
	return ret;
}

int scan_float (char const **ptr, float *x) {
	char *cur;
	errno = 0;
	*x = strtod (*ptr, &cur);
	int ret = !(cur == *ptr || errno) && isfinite (*x);
	*ptr = cur;
	return ret;
}

int scan2_unsigned (char const *s, unsigned *x) {
	return !scan_space (&s) && scan_unsigned (&s, x) && !*s;
}

int scan2_signed (char const *s, signed *x) {
	return !scan_space (&s) && scan_signed (&s, x) && !*s;
}

int scan2_float (char const *s, float *x) {
	return !scan_space (&s) && scan_float (&s, x) && !*s;
}

int scan_arr (
	char const **ptr, void *buf, unsigned n, unsigned size, scan_f *scan_one
) {
	for (unsigned i = 0;;) {
		if (!(*scan_one) (ptr, (char *) buf + i * size)) return 0;
		else if (++i == n) return 1;
		else if (!scan_space (ptr)) return 0;
	}
}

int scan_iarr (char const **ptr, signed buf[], unsigned n) {
	return scan_arr (ptr, buf, n, sizeof (signed), (scan_f *) &scan_signed);
}

int scan_uarr (char const **ptr, unsigned buf[], unsigned n) {
	return scan_arr (ptr, buf, n, sizeof (unsigned), (scan_f *) &scan_unsigned);
}

int scan_farr (char const **ptr, float buf[], unsigned n) {
	return scan_arr (ptr, buf, n, sizeof (float), (scan_f *) &scan_float);
}

unsigned more (unsigned a, unsigned b) {
	return a > b ? a : b;
}

int scan_block (
	struct cursor *curs, unsigned size, unsigned n0,
	unsigned *cnt, void **ret, scan_f *scan_one
) {
	void *buf = NULL, *tmp;
	unsigned n = n0, i = 0;

	if (n && !(buf = malloc (n * size)))
		{ cursor_err (curs, "E: failed malloc()"); goto err; }
	while (1) {
		if (!cursor_step (curs))
			{ cursor_err (curs, "E: premature EOF"); goto err; }
		else if (!curs->len) {
			if (cnt || i == n) goto retn;
			else { cursor_err (curs, "E: premature blank in block"); goto err; }
		} else if (!cnt && i == n)
			{ cursor_err (curs, "E: non-blank tail in block"); goto err; }
		else if (scan_space (CUR_PTR (curs)))
			{ cursor_err (curs, "E: leading whitespace"); goto err; }
		else while (1) {
			if (i == n) {
				if (cnt) n = SOME_MORE (n);
				else { cursor_err (curs, "E: oversized block"); goto err; }
				if ((tmp = realloc (buf, n * size))) buf = tmp;
				else { cursor_err (curs, "E: failed realloc()"); goto err; }
			}
			if ((*scan_one) (CUR_PTR (curs), (char *) buf + i * size)) ++i;
			else { cursor_err (curs, "E: failed to parse stream"); goto err; }
			if (!scan_space (CUR_PTR (curs))) {
				if (!*curs->ptr) break;
				else { cursor_err (curs, "E: malformed whitespace"); goto err; }
			} else if (!*curs->ptr)
				{ cursor_err (curs, "E: trailing whitespace"); goto err; }
		}
	}

	retn:
	if (cnt) {
		*cnt = i;
		if ((tmp = realloc (buf, i * size))) buf = tmp;
	}
	*ret = buf;
	return 1;
	err:
	if (buf) free (buf);
	return 0;
}

#ifdef BENCHMARK
unsigned ns_delta (struct timespec const *t1, struct timespec const *t2) {
	return (unsigned) (t2->tv_nsec - t1->tv_nsec) +
		1000000000 * (t2->tv_sec - t1->tv_sec);
}

void stat_step (float stat[2], float x) {
	stat[0] += x;
	stat[1] += x * x;
}

void stat_write (float const stat[2], unsigned const n[2]) {
	printf ("%.0f(%.0f)", stat[0] / (n[0] * n[1]), sqrt (
		(stat[1] - stat[0] * stat[0] / n[0]) / (n[0] - 1)
	) / n[1]);
}
#endif

