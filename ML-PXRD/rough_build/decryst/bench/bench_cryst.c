#define _POSIX_C_SOURCE 199309L  // clock_gettime(), etc.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rng.h"
#include "cryst.h"
#include "bench_cryst.h"
#include "bench_utils.h"

#define ALG_SAMPLE 100
#define ALG_TRIAL 100

typedef int work_f (crystal *);

struct algorithm {
	char const *desc;
	work_f *work;
};

static int xyzs_scan (float xyzs[][3], char const *const *argv, unsigned n) {
	for (unsigned i = 0; i < n; ++i) if (sscanf (
		argv[i], "%f,%f,%f", xyzs[i] + 0, xyzs[i] + 1, xyzs[i] + 2
	) != 3) return 0;
	return 1;
}

static int test_spec (
	crystal *cr, float xyzs[][3], struct algorithm *alg, unsigned n
) {
	cryst_set (cr, xyzs);
	if (!(*alg->work) (cr)) return 0;
	cryst_ack (cr, 1);
	printf ("%s: ", alg->desc);
	cryst_write (cr);
	return 1;
}

static int test_time (crystal *cr, rng_t *rng, struct algorithm *alg) {
	struct timespec t1, t2;
	float stat[2] = { 0.0 };
	unsigned n[2] = { ALG_SAMPLE, ALG_TRIAL };
	for (unsigned i = 0; i < ALG_SAMPLE; ++i) {
		if (clock_gettime (CLOCK_MONOTONIC, &t1)) return 0;
		for (unsigned j = 0; j < ALG_TRIAL; ++j) {
			cryst_step (cr, rng_float (rng));
			if (!(*alg->work) (cr)) return 0;
			cryst_ack (cr, 1);
		}
		if (clock_gettime (CLOCK_MONOTONIC, &t2)) return 0;
		stat_step (stat, ns_delta (&t1, &t2));
	}
	printf ("%s: ", alg->desc);
	stat_write (stat, n);
	printf ("\n");
	return 1;
}

int main (int argc, char const *const *argv) {
	struct algorithm aCr = { "cryst_eval", &cryst_eval },
		aBc = { "cryst_eval_bc", &cryst_eval_bc };
	rng_t *rng;
	struct brin_t *brin;
	crystal *cr;
	float (*xyzs)[3];
	unsigned nAtom;
	int ret = 0;

	if (!(rng = rng_mk2 ())) goto rng_err;
	if (!(brin = brin_read (stdin))) goto cr_err;
	nAtom = brin->crin->nAtom;
	cr = brin_eval (brin, rng);
	brin_fin (brin);
	if (!cr) goto cr_err;
	if (!(cryst_dof (cr) && (
		xyzs = malloc (nAtom * 3 * sizeof (float))
	))) goto xyzs_err;

	if (argc == 1) {
		if (!cryst_eval (cr)) goto retn_err;
		cryst_ack (cr, 1);
		ret = test_time (cr, rng, &aBc) && test_time (cr, rng, &aCr);
	} else {
		if (!(
			argc == nAtom + 1 && xyzs_scan (xyzs, argv + 1, nAtom)
		)) goto retn_err;
		ret = test_spec (cr, xyzs, &aBc, nAtom) &&
			test_spec (cr, xyzs, &aCr, nAtom);
	} retn_err:
	free (xyzs); xyzs_err:
	bryst_fin (cr); cr_err:
	free (rng); rng_err:
	return !ret;
}

