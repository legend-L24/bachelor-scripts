enum { CMD_FIN, CMD_INIT, CMD_COMP, CMD_DUMP, CMD_BEST, CMD_LOAD, CMD_CNT };

struct msg_init {
	uint8_t cmd[4];
	// (delta; tau, num of cores).
	uint32_t ff[1], uu[2];
};

struct msg_comp {
	uint8_t cmd[4];
	// (a, b, d, e, mag, r, latest s; idx of core).
	uint32_t ff[7], uu[1];
};

struct msg_cret {
	uint8_t cmd[4];
	// (sum, sq sum of deviation, final s; acceptance cnt).
	uint32_t ff[3], uu[1];
};

