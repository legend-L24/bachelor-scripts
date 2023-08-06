struct rbnode {
	// `cross' for (node, peer) <-> (peer, node) lookup.
	struct rbnode *parent, *left, *right, *cross;
	// `peer' is also the key.
	unsigned red, peer;
	signed score;
};

struct rbtree {
	struct rbnode nil, *root;
};

struct rbforest {
	struct rbtree *trees;
	// Total collision score for all unordered pairs.
	signed score;
	// cnt: number of unordered pairs.
	unsigned cnt, n;
};

// Returns NULL iff tree is empty.
extern struct rbnode *rbtree_first (struct rbtree const *tree);
// Returns NULL iff already maximum.
extern struct rbnode *rbtree_next
	(struct rbtree const *tree, struct rbnode *node);
// Returns NULL iff not found.
extern struct rbnode *rbtree_get (struct rbtree const *tree, unsigned x);

extern void rbtree_mk (struct rbtree *tree);
extern void rbtree_clear (struct rbtree *tree);
// Assumes no duplicate.
extern void rbtree_add (struct rbtree *tree, struct rbnode *node);
extern void rbtree_rm (struct rbtree *tree, struct rbnode *node);

// Initialises a forest of `n' rbtrees; returns 1 on success, 0 on failure.
extern int rbforest_mk (struct rbforest *forest, unsigned n);
extern void rbforest_fin (struct rbforest *forest);

// Adds unordered collision pair {x1, x2} (x1 != x2; x1, x2 < n) with `score'
// into `forest', assuming no duplicate; returns 1 on success, 0 on failure.
extern int rbforest_add
	(struct rbforest *forest, unsigned x1, unsigned x2, signed score);
// Removes all collision pairs containing `x' from `forest'.
extern void rbforest_rm (struct rbforest *forest, unsigned x);
// Writes sorted pairs.
extern void rbforest_pairs (struct rbforest const *forest, unsigned pairs[][2]);

