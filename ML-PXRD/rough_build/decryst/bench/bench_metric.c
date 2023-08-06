#define _POSIX_C_SOURCE 199309L  // clock_gettime(), etc.

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rng.h"
#include "metric.h"
#include "utils.h"
#include "bench_metric.h"
#include "bench_utils.h"

#define ID_TRIAL 1000000
#define ALG_SAMPLE 1000
#define ALG_TRIAL 1000
#define ALG_CNT (sizeof (Algs) / sizeof (Algs[0]))

typedef float work_f (metric const *, float const[3]);

struct algorithm {
	char const *desc;
	work_f *work;
};

struct parameter {
	float abc[3], angles[3], xyz[3];
	metric *met;
};

struct algorithm const Algs[] = {
	{ "metric_get_brute", &metric_get_brute },
	{ "metric_get_slow_bm", &metric_get_slow_bm },
	{ "metric_get_slow", &metric_get_slow },
	{ "metric_get_fast_bm", &metric_get_fast_bm },
	{ "metric_get_fast", &metric_get_fast }
};

// Requires all three pairs of faces of the generated cell to be overlapping
// under orthogonal projection, cf. Subsection 2.3.1 in the thesis for the
// motivation of this criterion.  Following is the proof.
//
// Let the direct and reciprocal bases be (a1, a2, a3) and (b1, b2, b3),
// respectively, then by the definition of a reciprocal basis we have
//   transpose((b1 b2 b3)) (a1 a2 a3) = transpose(B) A = I_3.
// So when using the direct basis for the vector space, the matrix
// representation for the reciprocal basis is
//   F = inverse(A) B = transpose(B) B [= inverse(transpose(A) A) = inverse(G)],
// where G and F are matrix representations of metric tensors for the direct and
// reciprocal spaces, respectively.
//
// Since the reciprocal basis contains the vectors normal to the three pairs of
// faces of the cell, the vectors are also direction vectors of the orthogonal
// projections in the proposition to be proved, and only need proper scaling in
// order to be the displacement vectors of the projections.  Noticing that the
// cell is a cube when using the direct basis, the displacement vectors can be
// computed from F:
//   ( f11 f21 f31 )      (    1    f21/f22 f31/f33 )
//   ( f12 f22 f32 )  ->  ( f12/f11    1    f32/f33 )
//   ( f13 f23 f33 )      ( f13/f11 f23/f22    1    )
//
// Noticing that the faces are squares when using the direct basis, the
// overlapping condition is satisfied iff the six latitudinal components of
// the three displacement vectors are in [-1, 1], so we only need to check
// whether the six elements in the above-mentioned matrix are all in said range.
// The formula for the elements of F can be looked up in the International
// Tables of Crystallography, and the rest is programming work.
//
// Q.E.D.

static int metric_chk (float const abc[3], float const angles[3]) {
	static unsigned const idx[3][3] = {{ 0, 1, 2 }, { 2, 0, 1 }, { 1, 2, 0 }};
	float cs[3], s2s[3];
	unsigned i, j;
	for (i = 0; i < 3; ++i) {
		cs[i] = cos (angles[i] * PI / 180.0);
		s2s[i] = 1 - cs[i] * cs[i];
	}
	for (i = 0; i < 3; ++i) {
		unsigned const *cur = idx[i];
		float r = cs[cur[0]] * cs[cur[1]] - cs[cur[2]];
		for (j = 0; j < 2; ++j) if (fabs (
			(r * abc[cur[j]]) / (abc[cur[!j]] * s2s[cur[j]])
		) > 1.0) return 0;
	}
	return 1;
}

static int param_mk (struct parameter *param, rng_t *rng) {
	float gap;
	unsigned i;
	do for (i = 0; i < 3; ++i) {
		param->abc[i] = 1.0 + rng_float (rng);
		param->angles[i] = 60.0 + 60.0 * rng_float (rng);
	} while (!metric_chk (param->abc, param->angles));
	do {
		gap = 0.0;
		for (i = 0; i < 3; ++i) {
			param->xyz[i] = 2.0 * rng_float (rng) - 1.0;
			gap += fabs (param->xyz[i] - round (param->xyz[i]));
		}
	} while (gap >= 0.5);
	return (param->met = metric_mk (param->abc, param->angles)) != NULL;
}

static int test_id (rng_t *rng) {
	struct parameter param;
	float ret[ALG_CNT];
	unsigned i, j;
	printf ("# Identity of proposed algorithms:\n");
	for (i = 0; i < ID_TRIAL; ++i) {
		if (!param_mk (&param, rng)) return 0;
		for (j = 0; j < ALG_CNT; ++j) ret[j] =
			(*Algs[j].work) (param.met, param.xyz);
		for (j = 1; j < ALG_CNT; ++j) if (fabs (ret[j] - ret[0]) > 1e-5) {
			printf (
				"Failure at %u\nabc: %.3f, %.3f, %.3f\n"
				"angles: %.1f, %.1f, %.1f\nxyz: %.3f, %.3f, %.3f\nret:",
				i, param.abc[0], param.abc[1], param.abc[2],
				param.angles[0], param.angles[1], param.angles[2],
				param.xyz[0], param.xyz[1], param.xyz[2]
			);
			for (j = 0; j < ALG_CNT; ++j) printf
				("%s%.3f", j ? ", " : " ", ret[j]);
			printf ("\n");
			break;
		}
		free (param.met);
	}
	return 1;
}

static int test_perf (rng_t *rng) {
	struct timespec t1, t2;
	struct parameter params[ALG_TRIAL];
	float stat[ALG_CNT][2] = {{ 0.0 }};
	unsigned n[2] = { ALG_SAMPLE, ALG_TRIAL }, perm[ALG_CNT], i, j, k;
	printf ("# Comparison of performance:\n");

	for (i = 0; i < ALG_CNT; ++i) perm[i] = i;
	for (i = 0; i < ALG_SAMPLE; ++i) {
		for (j = 0; j < ALG_TRIAL; ++j) {
			if (!param_mk (params + j, rng)) goto err;
		}
		rng_shuf (rng, perm, ALG_CNT);
		for (j = 0; j < ALG_CNT; ++j) {
			work_f *work = Algs[perm[j]].work;
			if (clock_gettime (CLOCK_MONOTONIC, &t1)) goto err;
			for (k = 0; k < ALG_TRIAL; ++k) (*work)
				(params[k].met, params[k].xyz);
			if (clock_gettime (CLOCK_MONOTONIC, &t2)) goto err;
			stat_step (stat[perm[j]], ns_delta (&t1, &t2));
		}
		for (j = 0; j < ALG_TRIAL; ++j) free (params[j].met);
	}
	for (i = 0; i < ALG_CNT; ++i) {
		printf ("%s: ", Algs[i].desc);
		stat_write (stat[i], n);
		printf ("\n");
	}

	return 1;
	err:
	for (i = 0; i < ALG_TRIAL; ++i) {
		if (params[i].met) free (params[i].met);
		else break;
	}
	return 0;
}

int main (void) {
	rng_t *rng;
	uint64_t seed;
	int ret;
	if (!(seed_mk (&seed) && (rng = rng_mk (seed)))) return 1;
	printf ("# Seed: 0x%016" PRIx64 "\n", seed);
	ret = !(test_id (rng) && test_perf (rng));
	free (rng);
	return ret;
}

