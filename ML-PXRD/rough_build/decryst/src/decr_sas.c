#define _XOPEN_SOURCE 500  // snprintf().

#include <math.h>
#include <pthread.h>
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
#include "int_netsa.h"

#define ADDR_MAX 32
#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
#define ERR_PROC(RET) \
	((void) (desc && fprintf (stderr, "E: %s in %s()\n", desc, __func__)), RET)

struct worker_t {
	struct crin_t const *crin;
	void *ctx, *buf;
	pthread_t thread;
	unsigned const (*sizes)[2];
	unsigned idx;
	int ret;
};

static int addr_print (char buf[], unsigned len, unsigned idx) {
	int rc = snprintf (buf, len, "inproc://%u", idx);
	return rc > 0 && rc < len;
}

static unsigned crin_dof (struct crin_t const *crin) {
	rng_t *rng;
	crystal *cr;
	char const *desc = NULL;
	unsigned dof = 0;
	if (!(rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
	if (!(cr = crin_eval (crin, rng))) ERR_GOTO ("failed crin_eval()", cr_err)
	if (!(dof = cryst_dof (cr))) fprintf
		(stderr, "E: dof == 0 in %s()\n", __func__);
	cryst_fin (cr); cr_err:
	free (rng); rng_err:
	return ERR_PROC(dof);
}

static int recv_chk (
	void *client, void *buf_, uint8_t *cmd_,
	unsigned const sizes[CMD_CNT + 1][2], unsigned n
) {
	int rc = zmq_recv (client, buf_, n * sizes[CMD_CNT][0], 0);
	if (rc < 1) return 0;
	uint8_t *buf = buf_, cmd = *buf;
	if (cmd >= CMD_CNT) return 0;
	unsigned size = sizes[cmd][0], i;
	if (rc != n * size) return 0;
	for (i = 1; i < n; ++i) if (*(buf + i * size) != cmd) return 0;
	*cmd_ = cmd;
	return 1;
}

static int send_chk (void *buf_, unsigned const size, unsigned n) {
	uint8_t const *buf = buf_;
	for (unsigned i = 0; i < n; ++i, buf += size) if (!buf[1]) return 0;
	return 1;
}

static int worker_mk (
	rng_t **rng, crystal **cr, sa_comp **comp, struct crin_t const *crin
) {
	char const *desc = NULL;
	if (!(*rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
	if (!(*cr = crin_eval (crin, *rng))) ERR_GOTO
		("failed crin_eval()", cr_err)
	if (!(*comp = sa_comp_mk (*cr, *rng))) ERR_GOTO
		("failed sa_comp_mk()", comp_err)
	if (!cryst_eval (*cr)) ERR_GOTO ("failed cryst_eval()", retn_err)
	cryst_ack (*cr, 1);
	return 1; retn_err:
	sa_comp_fin (*comp); comp_err:
	cryst_fin (*cr); cr_err:
	free (*rng); rng_err:
	return ERR_PROC(0);
}

// Seems to somehow trigger `-Wmaybe-uninitialized'?
// <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=79768>.
static void worker_fin (rng_t *rng, crystal *cr, sa_comp *comp) {
	sa_comp_fin (comp);
	cryst_fin (cr);
	free (rng);
}

static void *worker_main (struct worker_t *worker) {
	rng_t *rng;
	crystal *cr;
	sa_comp *comp;
	void *client;
	char const *desc = NULL;
	char addr[ADDR_MAX];
	float const *scores, *best;
	float ctl[8], stat[3];
	uint8_t *cb = worker->buf;
	uint32_t *ub = worker->buf;
	unsigned const (*sizes)[2] = worker->sizes, *bbuf;
	unsigned n[4], i;

	if (!addr_print (addr, sizeof (addr), worker->idx)) ERR_GOTO
		("failed snprintf()", client_err)
	if (!(client = zmq_socket (worker->ctx, ZMQ_REP))) ERR_GOTO
		("failed zmq_socket()", client_err)
	if (zmq_bind (client, addr)) ERR_GOTO ("failed zmq_bind()", worker_err)
	int flag = worker_mk (&rng, &cr, &comp, worker->crin);
	if (flag) {
		scores = cryst_scores (cr);
		best = sa_comp_best (comp);
		bbuf = sa_comp_bbuf (comp);
	}

	while (1) {
		if (zmq_recv (
			client, worker->buf, sizes[CMD_CNT][0], 0
		) == -1) ERR_GOTO ("failed zmq_recv()", retn_err)

		if (flag) switch (cb[0]) {
		case CMD_COMP:
			memcpy (ctl, ub + 1, 7 * sizeof (uint32_t));
			nltohf2 (ctl, 7);
			n[1] = ntohl (ub[8]);
			cb[1] = arefinite (ctl, 7) && ctl[6] >= 0.0 && n[1] < n[2] &&
				sa_comp_get (comp, n, ctl, n + 3, stat, 1);
			if (cb[1]) {
				memcpy (ub + 1, stat, 3 * sizeof (uint32_t));
				hftonl2 (ub + 1, 3);
				ub[4] = htonl (n[3]);
			}
			break;
		case CMD_DUMP:
			cb[1] = 1;
			for (i = 0; i < 3; ++i) ub[i + 1] = hftonl (scores[i]);
			cryst_dump (cr, ub + 4);
			break;
		case CMD_BEST:
			cb[1] = 1;
			for (i = 0; i < 3; ++i) ub[i + 1] = hftonl (best[i]);
			memcpy (ub + 4, bbuf, sizes[CMD_CNT][1]);
			break;
		case CMD_LOAD:
			if ((
				cb[1] = cryst_load (cr, ub + 1) && cryst_eval (cr)
			)) cryst_ack (cr, 1);
			break;
		case CMD_INIT:
			ctl[7] = nltohf (ub[1]);
			n[0] = ntohl (ub[2]);
			n[2] = ntohl (ub[3]);
			cb[1] = isfinite (ctl[7]) && ctl[7] > 1.0 && n[0] && n[2];
			break;
		default:
			break;
		} else cb[1] = 0;

		if (zmq_send (
			client, worker->buf, sizes[cb[0]][1], 0
		) == -1) ERR_GOTO ("failed zmq_send()", retn_err)
		if (cb[0] == CMD_FIN) {
			if (flag) break;
			else goto retn_err;
		}
	}

	worker->ret = 1; retn_err:
	if (flag) { worker_fin (rng, cr, comp); } worker_err:
	if (!mymq_close (client)) { worker->ret = 0; } client_err:
	return ERR_PROC(NULL);
}

static struct hosts_t *workers_run (
	void *ctx, struct crin_t const *crin, struct worker_t workers[],
	unsigned const sizes[CMD_CNT + 1][2], unsigned n
) {
	struct ihost_t *ihosts;
	struct hosts_t *hosts;
	char const *desc = NULL;
	char addr[ADDR_MAX];
	unsigned i;

	for (i = 0; i < n; ++i) {
		workers[i].crin = NULL;
		if (!(workers[i].buf = malloc (sizes[CMD_CNT][0]))) ERR_GOTO
			("failed malloc()", ihosts_err)
		workers[i].crin = crin;
		workers[i].ctx = ctx;
		workers[i].sizes = sizes;
		workers[i].idx = i;
		if (pthread_create (
			&workers[i].thread, NULL,
			(void *(*) (void *)) &worker_main, workers + i
		)) {
			free (workers[i].buf);
			workers[i].crin = NULL;
			ERR_GOTO ("failed pthread_create()", ihosts_err)
		}
	}
	if (!(ihosts = calloc (n, sizeof (struct ihost_t)))) ERR_GOTO
		("failed calloc()", ihosts_err)
	for (i = 0; i < n; ++i) {
		if (!addr_print (addr, sizeof (addr), i)) ERR_GOTO
			("failed snprintf()", retn_err)
		if (!(ihosts[i].socket = zmq_socket (ctx, ZMQ_REQ))) ERR_GOTO
			("failed zmq_socket()", retn_err)
		if (zmq_connect (ihosts[i].socket, addr)) ERR_GOTO
			("failed zmq_connect()", retn_err)
		ihosts[i].n = 1;
	}

	hosts = hosts_mk (ihosts, n);
	free (ihosts);
	return hosts; retn_err:
	for (i = 0; i < n; ++i) {
		if (ihosts[i].socket) mymq_close (ihosts[i].socket);
		else break;
	}
	free (ihosts); ihosts_err:
	return ERR_PROC(NULL);
}

static int workers_chk (struct worker_t workers[], unsigned n) {
	int ret = 1;
	for (unsigned i = 0; i < n; ++i) {
		if (!workers[i].crin) break;
		if (pthread_join (workers[i].thread, NULL) || !workers[i].ret) ret = 0;
		free (workers[i].buf);
	}
	return ret;
}

int main (int argc, char const *const *argv) {
	struct worker_t *workers;
	struct hosts_t *hosts;
	struct crin_t *crin;
	void *ctx, *client, *buf;
	char const *desc = NULL;
	uint8_t cmd;
	unsigned sizes[CMD_CNT + 1][2] = {
		{ 2 * sizeof (uint8_t), 0 },
		{ sizeof (struct msg_init), 2 * sizeof (uint8_t) },
		{ sizeof (struct msg_comp), sizeof (struct msg_cret) },
		{ 2 * sizeof (uint8_t), 0 },
		{ 2 * sizeof (uint8_t), 0 },
		{ 0, 2 * sizeof (uint8_t) },
		{ 0 }
	}, n;
	int flags[2] = { 0 };

	if (!(argc == 3 && scan2_unsigned (argv[1], &n))) {
		fprintf (stderr, "Usage: decr_sas num_of_cores zmq_addr < cryst.cr\n");
		goto workers_err;
	}
	if (!(workers = calloc (n, sizeof (struct worker_t)))) ERR_GOTO
		("failed calloc() for `workers'", workers_err)
	if (!(ctx = zmq_ctx_new ())) ERR_GOTO ("failed zmq_ctx_new()", ctx_err)
	if (!(client = zmq_socket (ctx, ZMQ_REP))) ERR_GOTO
		("failed zmq_socket()", client_err)
	if (zmq_bind (client, argv[2])) ERR_GOTO ("failed zmq_bind()", crin_err)
	if (!(crin = crin_read (stdin))) ERR_GOTO ("failed crin_read()", crin_err)
	if (!(sizes[CMD_CNT][1] = crin_dof (crin))) ERR_GOTO ("dof == 0", hosts_err)
	sizes[CMD_DUMP][1] = sizes[CMD_BEST][1] =
		(sizes[CMD_CNT][1] + 4) * sizeof (uint32_t);
	sizes[CMD_LOAD][0] = (sizes[CMD_CNT][1] + 1) * sizeof (uint32_t);
	sizes[CMD_CNT][0] = sizes[CMD_DUMP][1] < sizes[CMD_COMP][0] ?
		sizes[CMD_COMP][0] : sizes[CMD_DUMP][1];
	sizes[CMD_CNT][1] *= sizeof (uint32_t);
	if (!(hosts = workers_run (ctx, crin, workers, sizes, n))) goto hosts_err;
	if (!(buf = calloc (n, sizes[CMD_CNT][0]))) ERR_GOTO
		("failed calloc() for `buf'", buf_err)

	if (!recv_chk (client, buf, &cmd, sizes, n) || cmd != CMD_INIT) ERR_GOTO
		("error with initial request", retn_err)
	if (
		!vec_send (hosts, buf, sizes[cmd][0]) ||
		!vec_recv (hosts, buf, sizes[cmd][1]) ||
		zmq_send (client, buf, n * sizes[cmd][1], 0) == -1
	) ERR_GOTO ("error with initial reply", retn_err)
	crin_fin (crin); crin = NULL;
	while (1) {
		flags[1] = send_chk (buf, sizes[cmd][1], n);
		if (
			!recv_chk (client, buf, &cmd, sizes, n) ||
			cmd >= CMD_CNT || cmd == CMD_INIT ||
			!(flags[1] || cmd == CMD_FIN)
		) ERR_GOTO ("error with request", retn_err)
		if (
			!vec_send (hosts, buf, sizes[cmd][0]) ||
			!vec_recv (hosts, buf, sizes[cmd][1])
		) ERR_GOTO ("error with reply", retn_err)
		if (cmd == CMD_FIN) break;
		if (zmq_send (client, buf, n * sizes[cmd][1], 0) == -1) ERR_GOTO
			("error with reply", retn_err)
	}

	flags[0] = 1; retn_err:
	free (buf); buf_err:
	hosts_fin (hosts); hosts_err:
	if (crin) { crin_fin (crin); } crin_err:
	if (!mymq_close (client)) { flags[0] = 0; } client_err:
	if (zmq_ctx_destroy (ctx)) { flags[0] = 0; } ctx_err:
	if (!workers_chk (workers, n)) flags[0] = 0;
	free (workers); workers_err:
	return ERR_PROC(!flags[0]);
}

