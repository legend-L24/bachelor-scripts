#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zmq.h>
#include <arpa/inet.h>
#include "rng.h"
#include "cryst.h"
#include "optim.h"
#include "netio.h"
#include "utils.h"
#include "int_conf.h"
#include "int_netsa.h"

#define BUF_CNT 8
enum {
	BUF_PERM, BUF_BIAS, BUF_HEAD, BUF_INIT,
	BUF_COMP, BUF_CRET, BUF_DUMP, BUF_LOAD
};

struct main_res {
	struct hosts_t *hosts;
	rng_t *rng;
	sa_sched *sched;
	void *bufs[BUF_CNT], *ctx;
	float init[8];
	unsigned sizes[BUF_CNT + 1], n[4], last[2];
};

#define SCAN_VAL(SCAN, VAL) \
	if (!(++i < argc && SCAN (argv[i], VAL))) goto use_err;
#define CHK_DENY(DENY, NAME) \
	if (DENY) { fprintf (stderr, "E: invalid `%s'\n", NAME); ret = 0; }
static int init_scan (
	unsigned n[3], float init[5], char const *const *argv, int argc
) {
	unsigned i;
	int ret = 1;

	if (argc < 6) goto use_err;
	for (i = 1; i < argc; ++i) {
		if (!strcmp (argv[i], "-d")) {
			SCAN_VAL (scan2_float, init)
			CHK_DENY (init[0] <= 1.0, "delta")
		} else if (!strcmp (argv[i], "-k")) {
			SCAN_VAL (scan2_float, init + 1)
			CHK_DENY (init[1] <= 0.0, "kappa")
		} else if (!strcmp (argv[i], "-m")) {
			SCAN_VAL (scan2_unsigned, n + 2)
			CHK_DENY (!n[2], "mix")
		} else break;
	}

	if (!(i + 5 == argc &&
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
		"Usage: decr_sac [-d delta] [-k kappa] [-m mix] \\\n"
		"                tau lambda mean sd cryst.cr < hosts.txt\n"
	);
	return 0;
}

static void shuf_last (
	unsigned perm[], unsigned last[2], rng_t *rng, unsigned n, int init
) {
	unsigned i;
	if (init) for (i = 0; i < n; ++i) perm[i] = i;
	else last[1] = last[0];
	rng_shuf (rng, perm, n);
	for (i = 0;; ++i) if (perm[i] == n - 1) {
		last[0] = i;
		break;
	}
}

static int vec_chk (void const *buf_, unsigned size, unsigned n, uint8_t cmd) {
	uint8_t const *buf = buf_;
	for (unsigned i = 0; i < n; ++i, buf += size) {
		if (buf[0] != cmd || !buf[1]) return 0;
	}
	return 1;
}

static void head_fill (void *buf_, unsigned n, uint8_t cmd) {
	uint8_t (*buf)[2] = buf_;
	for (unsigned i = 0; i < n; ++i) buf[i][0] = cmd;
}

static void init_fill (
	struct msg_init inits[], unsigned const n[2], float delta
) {
	inits[0].cmd[0] = CMD_INIT;
	inits[0].ff[0] = hftonl (pow (delta, 1.0 / n[0]));
	inits[0].uu[0] = htonl (n[0]);
	inits[0].uu[1] = htonl (n[1]);
	for (unsigned i = 1; i < n[1]; ++i) memcpy
		(inits + i, inits + i - 1, sizeof (struct msg_init));
}

static void comp_fill_base (
	struct msg_comp comps[], unsigned const perm[], unsigned n
) {
	hftonl2 (comps[0].ff, 7);
	comps[0].uu[0] = htonl (perm[0]);
	for (unsigned i = 1; i < n; ++i) {
		memcpy (
			comps + i, comps + i - 1,
			sizeof (struct msg_comp) - sizeof (uint32_t)
		);
		comps[i].uu[0] = htonl (perm[i]);
	}
}

static int comp_fill0 (
	struct msg_comp comps[], sa_sched **sched,
	float const init[3], unsigned const perm[], unsigned n
) {
	if (!(*sched = sa_sched_mk (
		(float []) { ALPHA_TRIAL, BETA_TRIAL }, init, (float *) comps[0].ff
	))) return 0;
	comps[0].cmd[0] = CMD_COMP;
	comp_fill_base (comps, perm, n);
	return 1;
}

static int comp_fill (
	struct msg_comp comps[], sa_sched *sched,
	float const stat[4], unsigned const perm[], unsigned n
) {
	if (!sa_sched_step (sched, stat, (float *) comps[0].ff)) return 0;
	memcpy (comps[0].ff + 6, stat + 2, sizeof (float));
	comp_fill_base (comps, perm, n);
	return 1;
}

static int cret_proc (
	struct msg_cret const crets[],
	float stat[5], unsigned const n[2], unsigned last
) {
	float tmp[3];
	unsigned tau = n[0] * n[1], i;
	stat[4] = stat[0];
	memset (stat, 0, 4 * sizeof (float));
	for (i = 0; i < n[1]; ++i) {
		memcpy (tmp, crets[i].ff, sizeof (tmp));
		nltohf2 (tmp, 3);
		if (!arefinite (tmp, 3)) return 0;
		stat[0] += tmp[0];
		stat[1] += tmp[1];
		stat[3] += ntohl (crets[i].uu[0]);
	}
	stat[0] /= tau;
	stat[1] = sqrt (stat[1] / tau);
	stat[2] = nltohf (crets[last].ff[2]);
	stat[3] /= tau;
	return 1;
}

static int dump_proc (
	void const *dumps_, void *loads_, float bias[],
	rng_t *rng, unsigned const sizes[BUF_CNT + 1], unsigned n, float s
) {
	uint8_t const *dumps = dumps_ + sizeof (uint32_t);
	unsigned i;
	for (i = 0; i < n; ++i, dumps += sizes[BUF_DUMP]) {
		bias[i] = nltohf (*(uint32_t *) dumps);
		if (!isfinite (bias[i])) return 0;
	}

	float cap = bias[0];
	for (i = 1; i < n; ++i) if (bias[i] < cap) cap = bias[i];
	for (i = 0; i < n; ++i) bias[i] = exp ((cap - bias[i]) * s);
	wstoqs (bias, n);

	uint8_t *loads = loads_ + sizeof (uint32_t);
	dumps = dumps_ + 4 * sizeof (uint32_t);
	for (i = 0; i < n; ++i, loads += sizes[BUF_LOAD]) memcpy
		(loads, dumps + sizes[BUF_DUMP] * rng_bdice (rng, bias, n), sizes[BUF_CNT]);
	return 1;
}

static void dump_write (uint32_t const dumps[], unsigned const n[4]) {
	float score[2] = { nltohf (dumps[0]), 0.0 };
	unsigned size = n[3] + 4, best = 0, i, j;
	for (i = 1, j = size; i < n[1]; ++i, j += size) {
		if ((score[1] = dumps[j]) < score[0]) {
			best = i;
			score[0] = score[1];
		}
	}

	dumps += best * size;
	printf ("\n(%g %g %g)", nltohf (dumps[1]), nltohf (dumps[2]), nltohf (dumps[3]));
	dumps += 3;
	for (j = 0; j < n[3]; ++j) printf (" %g", nltohf (*++dumps));
	printf ("\n");
}

static void load_fill0 (void *loads_, unsigned size, unsigned n) {
	uint8_t *loads = loads_;
	for (unsigned i = 0; i < n; ++i, loads += size) loads[0] = CMD_LOAD;
}

#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
static int res_mk (struct main_res *res, char const *const *argv, int argc) {
	FILE *crf;
	crystal *cr;
	char const *desc = NULL;
	memcpy (res->init, (float []) {
		DELTA_TRIAL, KAPPA_TRIAL
	}, 2 * sizeof (float));
	memcpy (res->sizes, (unsigned []) {
		sizeof (unsigned), sizeof (float),
		2 * sizeof (uint8_t), sizeof (struct msg_init),
		sizeof (struct msg_comp), sizeof (struct msg_cret)
	}, 6 * sizeof (unsigned));
	res->n[2] = 0;

	if (!init_scan (res->n, res->init, argv, argc)) goto rng_err;
	if (!(res->rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
	if (!(crf = fopen (argv[argc - 1], "r"))) ERR_GOTO
		("failed fopen()", ctx_err)
	if (!(cr = cryst_read (res->rng, crf)))
		{ fclose (crf); ERR_GOTO ("failed cryst_read()", ctx_err) }
	if (fclose (crf))
		{ cryst_fin (cr); ERR_GOTO ("failed fclose()", ctx_err) }
	if (!(res->n[3] = cryst_dof (cr)))
		{ cryst_fin (cr); ERR_GOTO ("dof == 0", ctx_err) }
	cryst_fin (cr);
	if (!(res->ctx = zmq_ctx_new ())) ERR_GOTO
		("failed zmq_ctx_new()", ctx_err)
	if (!(res->hosts = hosts_read (res->ctx, stdin))) ERR_GOTO
		("failed hosts_read()", hosts_err)

	res->n[1] = res->hosts->n[1];
	if (!(res->n[2] || (
		res->n[2] = round (RATIO_TRIAL * res->n[1])
	))) res->n[2] = 1;
	res->sizes[BUF_DUMP] = (res->n[3] + 4) * sizeof (uint32_t);
	res->sizes[BUF_LOAD] = (res->n[3] + 1) * sizeof (uint32_t);
	res->sizes[BUF_CNT] = res->n[3] * sizeof (uint32_t);

	if (!bufs_mk (res->bufs, res->sizes, BUF_CNT, res->n[1])) ERR_GOTO
		("failed bufs_mk()", bufs_err)
	shuf_last (res->bufs[BUF_PERM], res->last, res->rng, res->n[1], 1);
	if (!comp_fill0 (
		res->bufs[BUF_COMP], &res->sched,
		res->init + 2, res->bufs[BUF_PERM], res->n[1]
	)) ERR_GOTO ("failed sa_sched_mk()", sched_err)
	init_fill (res->bufs[BUF_INIT], res->n, res->init[0]);
	load_fill0 (res->bufs[BUF_LOAD], res->sizes[BUF_LOAD], res->n[1]);

	return 1;
	free (res->sched); sched_err:
	bufs_fin (res->bufs, BUF_CNT); bufs_err:
	hosts_fin (res->hosts); hosts_err:
	zmq_ctx_destroy (res->ctx); ctx_err:
	free (res->rng); rng_err:
	if (desc) fprintf (stderr, "E: %s\n", desc);
	return 0;
}

static int res_fin (struct main_res *res) {
	int ret = 1;
	free (res->sched);
	bufs_fin (res->bufs, BUF_CNT);
	free (res->rng);
	if (!hosts_fin (res->hosts)) ret = 0;
	if (zmq_ctx_destroy (res->ctx)) ret = 0;
	return ret;
}

#define ERR_RETN(DESC) { desc = DESC; goto retn_err; }
#define DO_SEND(IDX) \
	if (!vec_send (res.hosts, res.bufs[IDX], res.sizes[IDX])) goto net_err;
#define DO_RECV(IDX, CMD) \
	if (!vec_recv (res.hosts, res.bufs[IDX], res.sizes[IDX])) goto net_err; \
	if (!vec_chk (res.bufs[IDX], res.sizes[IDX], res.n[1], CMD)) \
		ERR_RETN ("vec_chk()")
int main (int argc, char const *const *argv) {
	struct main_res res;
	char const *desc = "(wtf?)";
	int flags[2] = { 0 };

	if (!(res_mk (&res, argv, argc))) goto res_err;
	printf (
		"# decr_sa: (tau cores mix dof delta kappa lambda mean sd) ="
		" (%u %u %u %u %g %g %g %g %g)\n",
		res.n[0], res.n[1], res.n[2], res.n[3],
		res.init[0], res.init[1], res.init[2], res.init[3], res.init[4]
	);
	DO_SEND (BUF_INIT)
	DO_RECV (BUF_HEAD, CMD_INIT)
	printf ("[0]"); flags[1] = 1;

	for (unsigned i = 0, stable = 0;;) {
		DO_SEND (BUF_COMP)
		shuf_last (res.bufs[BUF_PERM], res.last, res.rng, res.n[1], 0);
		if (!(i % res.n[2]) && i) {
			dump_write (res.bufs[BUF_DUMP], res.n);
			printf ("[%u]", i);
		}
		DO_RECV (BUF_CRET, CMD_COMP)

		if (!cret_proc (
			res.bufs[BUF_CRET], res.init + 3, res.n, res.last[1]
		)) ERR_RETN ("cret_proc()")
		printf (" %g", res.init[3]);
		if (fabs (res.init[7] - res.init[3]) > res.init[1]) stable = 0;
		else if (++stable == STABLE_LIMIT) break;
		if (!comp_fill (
			res.bufs[BUF_COMP], res.sched,
			res.init + 3, res.bufs[BUF_PERM], res.n[1]
		)) ERR_RETN ("comp_fill()")

		if (!(++i % res.n[2])) {
			head_fill (res.bufs[BUF_HEAD], res.n[1], CMD_DUMP);
			DO_SEND (BUF_HEAD)
			DO_RECV (BUF_DUMP, CMD_DUMP)
			if (!dump_proc (
				res.bufs[BUF_DUMP], res.bufs[BUF_LOAD], res.bufs[BUF_BIAS],
				res.rng, res.sizes, res.n[1], res.init[5]
			)) ERR_RETN ("dump_proc()")
			DO_SEND (BUF_LOAD)
			DO_RECV (BUF_HEAD, CMD_LOAD)
		}
	}

	head_fill (res.bufs[BUF_HEAD], res.n[1], CMD_BEST);
	DO_SEND (BUF_HEAD)
	DO_RECV (BUF_DUMP, CMD_BEST)
	dump_write (res.bufs[BUF_DUMP], res.n);
	flags[1] = 0; flags[0] = 1; retn_err:
	if (flags[1]) {
		printf ("\n"); flags[1] = 0;
		fprintf (stderr, "E: failed %s\n", desc);
	}

	head_fill (res.bufs[BUF_HEAD], res.n[1], CMD_FIN);
	if (!vec_send (
		res.hosts, res.bufs[BUF_HEAD], res.sizes[BUF_HEAD]
	)) {
		flags[0] = 0;
		fprintf (stderr, "E: error with message transmission\n");
	} net_err:
	if (flags[1]) {
		printf ("\n"); flags[1] = 0;
		fprintf (stderr, "E: error with message transmission\n");
	}
	if (!res_fin (&res)) { flags[0] = 0; } res_err:
	return !flags[0];
}

