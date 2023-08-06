#include <stdint.h>
#include <stdio.h>

#define PI 3.14159265358979323846
// Maximum line width, including the line terminator and a `\0'.
#define LEN_MAX 1024

struct cursor {
	FILE *file;
	// Contents of current line, and the cursor.
	char buf[LEN_MAX], *ptr;
	// (index, length) of current line.
	unsigned idx, len;
};

typedef int scan_f (char const **, void *);

#define CUR_PTR(CURS) ((char const **) &CURS->ptr)
#define SOME_MORE(N) more (1.618 * (N), (N) + 1)
extern unsigned more (unsigned a, unsigned b);

// Float {,de}serialisation, also supporting +-nan, +-inf, -0 and denormals.
// The payload of (+-)nan is not guaranteed to be preserved.
extern uint32_t hftol (float x);
extern float hltof (uint32_t x);
extern uint32_t hftonl (float x);
extern float nltohf (uint32_t x);
// The above applied on a 4-byte aligned `buf'fer.
extern void hftonl2 (void *buf, unsigned n);
extern void nltohf2 (void *buf, unsigned n);
// Whether all numbers in `buf' are finite.
extern int arefinite (float const buf[], unsigned n);

// Initialises the cursor at line -1.
extern void cursor_mk (struct cursor *curs, FILE *file);
// Writes cursor position and error `desc'ription to stderr.
extern void cursor_err (struct cursor *curs, char const desc[]);
// Reads the next line, ignoring all lines beginning with `#'; returns 0 on EOF,
// or 1 otherwise.
extern int cursor_step (struct cursor *curs);

// Scans past consecutive whitespaces (' ' and '\t') beginning from `*ptr';
// returns 1 if `*ptr' points to whitespace, or 0 otherwise.
extern int scan_space (char const **ptr);
// All of these three functions return 1 on success, 0 on failure.
extern int scan_unsigned (char const **ptr, unsigned *x);
extern int scan_signed (char const **ptr, signed *x);
extern int scan_float (char const **ptr, float *x);
// Scan from a `s'tring containing only the number, forbidding any leading or
// trailing whitespace; return 1 on success, 0 on failure.
extern int scan2_unsigned (char const *s, unsigned *x);
extern int scan2_signed (char const *s, signed *x);
extern int scan2_float (char const *s, float *x);

// Scans a `buf'fer of length [size][n] (both `size' and `n' shall be > 0) with
// whitespace-delimited elements by applying `scan_one' on each element; returns
// 1 on success, 0 on failure.
extern int scan_arr
	(char const **ptr, void *buf, unsigned n, unsigned size, scan_f *scan_one);
// `scan_arr' instances with actual types; return 1 on success, 0 on failure.
extern int scan_iarr (char const **ptr, signed buf[], unsigned n);
extern int scan_uarr (char const **ptr, unsigned buf[], unsigned n);
extern int scan_farr (char const **ptr, float buf[], unsigned n);

// Reads a block of elements delimited by whitespaces and/or newlines, with the
// end marked by a blank line, disallowing leading or trailing whitespace on any
// line.  If `cnt' is NULL, the element count shall be exactly `n0'; otherwise
// memory will be pre-allocated for `n0' elements and grow as needed, and `*cnt'
// will be set to the final count on successful +return.  Returns 1 on success, 0
// on failure.
extern int scan_block (
	struct cursor *curs, unsigned size, unsigned n0,
	unsigned *cnt, void **ret, scan_f *scan_one
);

