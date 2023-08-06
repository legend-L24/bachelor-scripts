#include <stdint.h>
#include <stdlib.h>
#include <zmq.h>
#include "utils.h"
#include "netio.h"

typedef int mymq_sr_f (void *, void *, size_t, int);

int mymq_close (void *socket) {
	int linger = 100,
		ret = !zmq_setsockopt (socket, ZMQ_LINGER, &linger, sizeof (linger));
	return !zmq_close (socket) && ret;
}

static int ihosts_close (struct ihost_t const ihosts[], unsigned n) {
	unsigned i;
	int ret = 1;
	for (i = 0; i < n; ++i) if (!mymq_close (ihosts[i].socket)) ret = 0;
	return ret;
}

struct hosts_t *hosts_mk (struct ihost_t const ihosts[], unsigned n) {
	struct hosts_t *hosts;
	unsigned i;
	if (!(hosts = malloc (sizeof (struct hosts_t)))) goto hosts_err;
	if (!(hosts->items = calloc (n, sizeof (zmq_pollitem_t)))) goto items_err;
	if (!(hosts->hs = calloc (n, sizeof (unsigned)))) goto hs_err;

	for (hosts->n[1] = i = 0, hosts->n[0] = n; i < n; ++i) {
		hosts->items[i].socket = ihosts[i].socket;
		hosts->hs[i] = ihosts[i].n;
		hosts->n[1] += ihosts[i].n;
	}
	return hosts;
	free (hosts->hs); hs_err:
	free (hosts->items); items_err:
	free (hosts); hosts_err:
	ihosts_close (ihosts, n);
	return NULL;
}

int hosts_fin (struct hosts_t *hosts) {
	int ret = 1;
	for (unsigned i = 0; i < hosts->n[0]; ++i) {
		if (!mymq_close (hosts->items[i].socket)) ret = 0;
	}
	free (hosts->hs);
	free (hosts->items);
	free (hosts);
	return ret;
}

static int ihost_read (void *ctx, struct ihost_t *ihost, struct cursor *curs) {
	if (!(ihost->socket = zmq_socket (ctx, ZMQ_REQ)))
		{ cursor_err (curs, "E: failed zmq_socket()"); goto socket_err; }
	else if (scan_space (CUR_PTR (curs)))
		{ cursor_err (curs, "E: leading whitespace"); goto retn_err; }
	else if (!(scan_unsigned (CUR_PTR (curs), &ihost->n) && ihost->n))
		{ cursor_err (curs, "E: malformed `n'"); goto retn_err; }
	else if (
		!scan_space (CUR_PTR (curs)) || zmq_connect (ihost->socket, curs->ptr)
	) { cursor_err (curs, "E: failed zmq_connect()"); goto retn_err; }
	return 1; retn_err:
	mymq_close (ihost->socket); socket_err:
	return 0;
}

struct hosts_t *hosts_read (void *ctx, FILE *file) {
	struct cursor curs;
	struct ihost_t *ihosts = NULL, *tmp;
	struct hosts_t *hosts = NULL;
	unsigned n = 0, i = 0;
	int iret = 0;

	cursor_mk (&curs, file);
	for (; cursor_step (&curs); ++i) {
		if (i == n) {
			n = SOME_MORE (n);
			if ((
				tmp = realloc (ihosts, n * sizeof (struct ihost_t))
			)) ihosts = tmp;
			else { cursor_err (&curs, "E: failed realloc()"); goto err; }
		}
		if (!ihost_read (ctx, ihosts + i, &curs)) goto err;
	}

	iret = 1; err:
	if (ihosts) {
		if (iret) hosts = hosts_mk (ihosts, i);
		else ihosts_close (ihosts, i);
		free (ihosts);
	}
	return hosts;
}

static int vec_sendrecv (
	struct hosts_t *hosts, void *buf, unsigned size, short ev, mymq_sr_f *op
) {
	zmq_pollitem_t *items = hosts->items;
	unsigned *hs = hosts->hs, n = hosts->n[0], m = n, i;
	for (i = 0; i < n; ++i) items[i].events = ev;
	while (m) {
		unsigned idx = 0;
		int cnt = zmq_poll (items, n, -1);
		if (cnt == -1) return 0;
		for (i = 0, m -= cnt; i < n && cnt; idx += hs[i++]) {
			if (items[i].revents) {
				unsigned sz = hs[i] * size;
				int rc = (*op)
					(items[i].socket, (uint8_t *) buf + idx * size, sz, 0);
				if (rc == -1 || rc != sz) return 0;
				items[i].events = 0;
				--cnt;
			}
		}
	}
	return 1;
}

int vec_send (struct hosts_t *hosts, void const *buf, unsigned size) {
	return vec_sendrecv
		(hosts, (void *) buf, size, ZMQ_POLLOUT, (mymq_sr_f *) &zmq_send);
}

int vec_recv (struct hosts_t *hosts, void *buf, unsigned size) {
	return vec_sendrecv (hosts, buf, size, ZMQ_POLLIN, &zmq_recv);
}

int bufs_mk (void *bufs[], unsigned const sizes[], unsigned nb, unsigned bc) {
	unsigned i;
	for (i = 0; i < nb; ++i) if (!(
		bufs[i] = calloc (bc, sizes[i])
	)) goto err;
	return 1;
	err:
	for (i = 0; i < nb; ++i) {
		if (bufs[i]) free (bufs[i]);
		else break;
	}
	return 0;
}

void bufs_fin (void *bufs[], unsigned nb) {
	for (unsigned i = 0; i < nb; ++i) free (bufs[i]);
}

