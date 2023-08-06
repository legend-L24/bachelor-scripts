#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"
#include "cryst.h" // I use Val_spec from it
#include "optim.h"
#include "utils.h"
#include "int_conf.h"

#define SCAN_VAL(SCAN, VAL) \
	if (!(++i < argc && SCAN (argv[i], VAL))) goto use_err;
#define CHK_DENY(DENY, NAME) \
	if (DENY) { fprintf (stderr, "E: invalid `%s'\n", NAME); ret = 0; }
static int init_scan (
	unsigned n[1], float init[5], char const *const *argv, int argc
) {
	unsigned i;
	int ret = 1;

	if (argc < 5) goto use_err;
	for (i = 1; i < argc; ++i) {
		if (!strcmp (argv[i], "-d")) {
			SCAN_VAL (scan2_float, init)
			CHK_DENY (init[0] <= 1.0, "delta")
		} else if (!strcmp (argv[i], "-k")) {
			SCAN_VAL (scan2_float, init + 1)
			CHK_DENY (init[1] <= 0.0, "kappa")
		} else break;
	}

	if (!(
		i + 4 == argc &&
		scan2_unsigned (argv[i++], n) &&
		scan2_float (argv[i++], init + 2) &&
		scan2_float (argv[i++], init + 3) &&
		scan2_float (argv[i++], init + 4)
	)) goto use_err;
	CHK_DENY (!n[0], "tau")
	CHK_DENY (init[2] <= 0.0, "lambda")
	CHK_DENY (init[3] <= 0.0, "mean")
	CHK_DENY (init[4] <= 0.0, "sd")

	return ret; use_err:
	fprintf (
		stderr,
		"Usage: decr_lsa [-d delta] [-k kappa] tau lambda mean sd < cryst.cr\n"
	);
	return 0;
}


static void dump_write (float const *scores, uint32_t const buf[], unsigned n) {
	printf ("\n(%g %g %g)", scores[0], scores[1], scores[2]);
	for (unsigned i = 0; i < n; ++i) printf (" %g", nltohf (buf[i]));
	printf ("\n");
}
#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
int main (int argc, char const *const *argv) {
	rng_t *rng;
	crystal *cr;
	sa_sched *sched;
	sa_comp *comp;
	char const *desc = NULL;
	float const *scores;
	float init[8], ctl[8];
	uint32_t *buf;
	unsigned n[5];
	int flags[2] = { 0 };
	memcpy (init, (float []) { DELTA_TRIAL, KAPPA_TRIAL }, 2 * sizeof (float));
	memcpy (n + 1, (unsigned []) { 0, 1 }, 2 * sizeof (unsigned));

	if (!init_scan (n, init, argv, argc)) goto rng_err;
	if (!(rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
	if (!(cr =  cryst_read(rng, stdin))) ERR_GOTO
		("failed cryst_read()", cr_err)
	if (!cryst_eval (cr)) ERR_GOTO ("failed cryst_eval()", sched_err)
	if (!(n[3] = cryst_dof (cr))) ERR_GOTO ("dof == 0", sched_err) // I change n[3] to nPos[0]
	cryst_ack (cr, 1);
	scores = cryst_scores (cr);
	if (!(sched = sa_sched_mk (
		(float []) { ALPHA_TRIAL, BETA_TRIAL }, init + 2, ctl
	))) ERR_GOTO ("failed sa_sched_mk()", sched_err)
	ctl[7] = pow (init[0], 1.0 / n[0]);
	if (!(comp = sa_comp_mk (cr, rng))) ERR_GOTO
		("failed sa_somp_mk()", comp_err)
	if (!(buf = malloc (n[3] * sizeof (uint32_t)))) ERR_GOTO
		("failed malloc()", buf_err)
        /*
	printf (
		"# decr_lsa: (tau dof delta kappa lambda mean sd) ="
		" (%u %u %g %g %g %g %g): for cell parameter\n[0]",
		n[0], n[3], init[0], init[1], init[2], init[3], init[4]
	); flags[1] = 1;
// added to update cell parameter
	for (unsigned i = 0, stable = 0;;) {
		init[7] = init[3];
		if (!sa_comp_get (comp, n, ctl, n + 4, init + 3,0)) ERR_GOTO
			("failed sa_comp_get()", retn_err)
		init[3] /= n[0];
		init[4] = sqrt (init[4] / n[0]);
		init[6] = (float) n[4] / n[0];
		printf (" %g", init[3]);

		if (fabs (init[7] - init[3]) > init[1]) stable = 0;
		else if (++stable == STABLE_LIMIT) break;
		if (!sa_sched_step (sched, init + 3, ctl)) ERR_GOTO
			("failed sa_sched_step()", retn_err)
		if (!(++i % LINE_LIMIT)) {
			cryst_dump (cr, buf);
			dump_write (scores, buf, n[3]); 
			printf ("[%u]", i);
		}
		ctl[6] = init[5];
	}
	*/
	printf (
		"# decr_lsa: (tau dof delta kappa lambda mean sd) ="
		" (%u %u %g %g %g %g %g) for coordinate: \n[0]",
		n[0], n[3], init[0], init[1], init[2], init[3], init[4]
	); flags[1] = 1;

	for (unsigned i = 0, stable = 0;;) {
		init[7] = init[3];
		if (!sa_comp_get (comp, n, ctl, n + 4, init + 3,1)) ERR_GOTO
			("failed sa_comp_get()", retn_err)
		init[3] /= n[0];
		init[4] = sqrt (init[4] / n[0]);
		init[6] = (float) n[4] / n[0];
		printf (" %g", init[3]);

		if (fabs (init[7] - init[3]) > init[1]) stable = 0;
		else if (++stable == STABLE_LIMIT) break;
		if (!sa_sched_step (sched, init + 3, ctl)) ERR_GOTO
			("failed sa_sched_step()", retn_err)
		if (!(++i % LINE_LIMIT)) {
			cryst_dump (cr, buf);
			dump_write (scores, buf, n[3]);
			printf ("[%u]", i);
		}
		ctl[6] = init[5];
	}

	dump_write (sa_comp_best (comp), sa_comp_bbuf (comp), n[3]);
	flags[1] = 0; flags[0] = 1; retn_err:
	if (flags[1]) printf ("\n");
	free (buf); buf_err:
	sa_comp_fin (comp); comp_err:
	free (sched); sched_err:
	cryst_fin (cr); cr_err:
	free (rng); rng_err:
	if (desc) fprintf (stderr, "E: %s\n", desc);
	return !flags[0];
}
