#define _POSIX_C_SOURCE 199309L  // clock_gettime(), etc.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rng.h"
#include "cryst.h"
#include "bench_cryst.h"
#include "bench_utils.h"

#define ALG_COUNT 10000
#define ALG_SAMPLE 250
#define ALG_TRIAL 250

static int narrow_one (crystal const *cr, unsigned i, unsigned j, signed *ret) {
	return 1;
}

static int worker (crystal *cr, rng_t *rng, broad_f *broad, narrow_f *narrow) {
	cryst_step (cr, rng_float (rng));
	if (!cryst_eval_bb (cr, broad, narrow)) return 0;
	cryst_ack (cr, 1);
	return 1;
}

static int count_narrow (crystal *cr, rng_t *rng, narrow_f *narrow) {
	float const *scores = cryst_scores (cr);
	float stat[2] = { 0.0 };
	unsigned n[2] = { ALG_COUNT, 1 };
	for (unsigned i = 0; i < ALG_COUNT; ++i) {
		if (!worker (cr, rng, &bump_eval_bb, narrow)) return 0;
		stat_step (stat, scores[0]);
	}
	printf ("\t");
	stat_write (stat, n);
	fflush (stdout);
	return 1;
}

static int time_bump (
	crystal *cr, rng_t *rng, broad_f *broad, narrow_f *narrow
) {
	struct timespec t1, t2;
	float stat[2] = { 0.0 };
	unsigned n[2] = { ALG_SAMPLE, ALG_TRIAL };
	for (unsigned i = 0; i < ALG_SAMPLE; ++i) {
		if (clock_gettime (CLOCK_MONOTONIC, &t1)) return 0;
		for (unsigned j = 0; j < ALG_TRIAL; ++j) {
			if (!worker (cr, rng, broad, narrow)) return 0;
		}
		if (clock_gettime (CLOCK_MONOTONIC, &t2)) return 0;
		stat_step (stat, ns_delta (&t1, &t2));
	}
	printf ("\t");
	stat_write (stat, n);
	fflush (stdout);
	return 1;
}

int main (void) {
	rng_t *rng;
	struct brin_t *brin;
	crystal *cr;
	unsigned nAtom[2];
	int ret = 1;

	if (!(rng = rng_mk2 ())) goto rng_err;
	if (!(brin = brin_read (stdin))) goto cr_err;
	nAtom[0] = brin->crin->nAtom;
	cr = brin_eval (brin, rng);
	brin_fin (brin);
	if (!cr) goto cr_err;
	if (!(cryst_eval (cr))) goto retn_err;

	nAtom[1] = cryst_nball (cr);
	cryst_ack (cr, 1);
	printf (
		"\t%u\t%u\t%u\t%u", nAtom[1], nAtom[0], nAtom[1] * (nAtom[1] - 1) / 2,
		// nAtom * (nBall - nAtom) + nAtom * (nAtom - 1) / 2.
		nAtom[0] * nAtom[1] - nAtom[0] * (nAtom[0] + 1) / 2
	);
	ret = !(
		count_narrow (cr, rng, &narrow_one) &&
		count_narrow (cr, rng, &narrow_cryst) &&
		time_bump (cr, rng, NULL, NULL) &&
		time_bump (cr, rng, &bump_eval_bc, &narrow_one) &&
		time_bump (cr, rng, &bump_eval_bc, &narrow_cryst) &&
		time_bump (cr, rng, &bump_eval_bb, &narrow_one) &&
		time_bump (cr, rng, &bump_eval_bb, &narrow_cryst) &&
		time_bump (cr, rng, &bump_eval_orig, &narrow_orig)
	);

	retn_err:
	printf ("\n");
	bryst_fin (cr); cr_err:
	free (rng); rng_err:
	return ret;
}

