#include <stdio.h>
#include <zmq.h>

struct ihost_t {
	void *socket;
	unsigned n;  // Shall be > 0.
};

struct hosts_t {
	zmq_pollitem_t *items;
	// (number of cores on each host), (number of hosts, number of cores).
	// n[0] shall be > 0.
	unsigned *hs, n[2];
};

// Sets low LINGER on a `socket' and closes it; returns 1 on success, 0 on
// failure.
int mymq_close (void *socket);
// Returns NULL iff failed; on failure, sockets in `ihosts' will be closed.
extern struct hosts_t *hosts_mk (struct ihost_t const ihosts[], unsigned n);
// Closes all sockets in `hosts' and free all memory allocated with `hosts';
// returns 1 on success, 0 if any socket failed to close.
extern int hosts_fin (struct hosts_t *hosts);
// Creates REQ sockets, and connects them to REP endpoints in zmq `ctx', with
// the list of endpoints read from open `file'; returns NULL iff failed
// (including when nHost would be 0).
extern struct hosts_t *hosts_read (void *ctx, FILE *file);

// Sends/receives bytes in `buf'fer of length [nCore][size] (`nCore' and `size'
// shall be > 0); returns 1 on success, 0 on failure.
extern int vec_send (struct hosts_t *hosts, void const *buf, unsigned size);
extern int vec_recv (struct hosts_t *hosts, void *buf, unsigned size);

// Makes `nb' buffers, each of `bc' blocks, with block size `sizes' respectively
// for each buffer; returns 1 on success, 0 on failure.
extern int bufs_mk
	(void *bufs[], unsigned const sizes[], unsigned nb, unsigned bc);
extern void bufs_fin (void *bufs[], unsigned nb);

