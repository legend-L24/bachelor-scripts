#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
#include "cryst.h"
#include "optim.h"
#include "utils.h"

#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
int main (int argc, char const *const *argv) {
	rng_t *rng;
	crystal *cr;
	char const *desc = NULL;
	unsigned n;
	int ret = 1;

	if (!(argc == 2 && scan2_unsigned (argv[1], &n))) {
		fprintf (stderr, "Usage: decr_mcs num_of_samples < cryst.cr\n");
		goto rng_err;
	}
	if (!(rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
	if (!(cr = cryst_read (rng, stdin))) ERR_GOTO ("failed cryst_read()", cr_err)
        if (!peak_fun (cr)) ERR_GOTO ("failed peak_fun()", retn_err) // I merge two mistakes
	if (!(cryst_eval (cr))) ERR_GOTO ("failed cryst_eval()", retn_err)
	cryst_ack (cr, 1);
	if (!(mc_grid(cr))) ERR_GOTO ("failed mc_grid()", retn_err);
	/*
	if (cryst_dof (cr)) {
		if (n < 2) ERR_GOTO ("n < 2 with dof != 0", retn_err)
		if (!mc_stat (cr, rng, n, stat)) ERR_GOTO ("failed mc_stat()", retn_err)
		printf ("%g %g %g %g %g\n", stat[0], stat[1], stat[2], stat[3], stat[4]);
	} else {
		if (n) ERR_GOTO ("n != 0 with dof == 0", retn_err)
		float const *scores = cryst_scores (cr);
		printf ("- - %g %g %g\n", scores[0], scores[1], scores[2]);
	}
	*/
	ret = 0; retn_err:
	cryst_fin (cr); cr_err:
	free (rng); rng_err:
	if (desc) fprintf (stderr, "E: %s\n", desc);
	return ret;
}

