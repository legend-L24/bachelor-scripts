#include <stdlib.h>
#include "rbtree.h"

// Implementation of the rbtree modeled from CLRS (3rd edition).

static struct rbnode *rbnode_first (
	struct rbnode const *nil, struct rbnode *node
) {
	while (node->left != nil) node = node->left;
	return node;
}

struct rbnode *rbtree_first (struct rbtree const *tree) {
	struct rbnode const *nil = &tree->nil;
	return tree->root == nil ? NULL : rbnode_first (nil, tree->root);
}

struct rbnode *rbtree_next (struct rbtree const *tree, struct rbnode *node) {
	struct rbnode const *nil = &tree->nil;
	if (node->right != nil) return rbnode_first (nil, node->right);
	while (1) {
		if (node->parent == nil) return NULL;
		if (node == node->parent->left) return node->parent;
		node = node->parent;
	}
}

struct rbnode *rbtree_get (struct rbtree const *tree, unsigned x) {
	struct rbnode const *nil = &tree->nil;
	struct rbnode *node = tree->root;
	while (node != nil) {
		if (node->peer == x) return node;
		else node = x < node->peer ? node->left : node->right;
	}
	return NULL;
}

void rbtree_mk (struct rbtree *tree) {
	tree->nil = (struct rbnode) { .parent = &tree->nil, .red = 0 };
	tree->root = &tree->nil;
}

static void rbnode_clear (struct rbnode const *nil, struct rbnode *node) {
	if (node == nil) return;
	rbnode_clear (nil, node->left);
	rbnode_clear (nil, node->right);
	free (node);
}

void rbtree_clear (struct rbtree *tree) {
	rbnode_clear (&tree->nil, tree->root);
	tree->root = &tree->nil;
}

static void rbtree_graft (
	struct rbtree *tree, struct rbnode *node, struct rbnode *bro
) {
	if ((bro->parent = node->parent) == &tree->nil) tree->root = bro;
	else if (node == node->parent->left) node->parent->left = bro;
	else node->parent->right = bro;
}

static void rbtree_rotate (struct rbtree *tree, struct rbnode *node, int left) {
	struct rbnode *bro;
	if (left) {
		bro = node->right;
		(node->right = bro->left)->parent = node;
		bro->left = node;
	} else {
		bro = node->left;
		(node->left = bro->right)->parent = node;
		bro->right = node;
	}
	rbtree_graft (tree, node, bro);
	node->parent = bro;
}

void rbtree_add (struct rbtree *tree, struct rbnode *node) {
	struct rbnode *nil = &tree->nil, *parent = nil, *tmp = tree->root, *grand;
	int left = 0;
	while (tmp != nil) {
		parent = tmp;
		left = node->peer < parent->peer;
		tmp = left ? parent->left : parent->right;
	}
	node->red = 1;
	node->left = node->right = nil;
	if ((node->parent = parent) == nil) tree->root = node;
	else if (left) parent->left = node;
	else parent->right = node;

	while (parent->red) {
		left = parent == (grand = parent->parent)->left;
		tmp = left ? grand->right : grand->left;
		if (tmp->red) {
			parent->red = tmp->red = 0; grand->red = 1;
			parent = (node = grand)->parent;
		} else {
			if (node == (left ? parent->right : parent->left)) {
				rbtree_rotate (tree, parent, left);
				tmp = parent; parent = node; node = tmp;
			}
			parent->red = 0; grand->red = 1;
			rbtree_rotate (tree, grand, !left);
		}
	}
	tree->root->red = 0;
}

void rbtree_rm (struct rbtree *tree, struct rbnode *node) {
	struct rbnode *nil = &tree->nil, *orig = node, *parent, *tmp;
	int rl;
	if (node->left == nil) {
		rbtree_graft (tree, node, node->right);
		rl = node->red;
		node = node->right;
	} else if (node->right == nil) {
		rbtree_graft (tree, node, node->left);
		rl = node->red;
		node = node->left;
	} else {
		struct rbnode *broken = (tmp = rbnode_first (nil, node->right))->right;
		rl = tmp->red; tmp->red = node->red;
		(tmp->left = node->left)->parent = tmp;
		if ((parent = tmp->parent) == node) broken->parent = tmp;
		else {
			(parent->left = broken)->parent = parent;
			(tmp->right = node->right)->parent = tmp;
		}
		rbtree_graft (tree, node, tmp);
		node = broken;
	}
	if (rl) goto retn;

	while (!((parent = node->parent) == nil || node->red)) {
		rl = node == parent->left;
		tmp = rl ? parent->right : parent->left;
		if (tmp->red) {
			tmp->red = 0; parent->red = 1;
			tmp = rl ? tmp->left : tmp->right;
			rbtree_rotate (tree, parent, rl);
		}
		if (!(tmp->left->red || tmp->right->red)) {
			tmp->red = 1;
			node = parent;
		} else {
			if (!(rl ? tmp->right : tmp->left)->red) {
				struct rbnode *child = rl ? tmp->left : tmp->right;
				child->red = 0; tmp->red = 1;
				rbtree_rotate (tree, tmp, !rl);
				tmp = child;
			}
			tmp->red = parent->red;
			(rl ? tmp->right : tmp->left)->red = parent->red = 0;
			rbtree_rotate (tree, parent, rl);
			node = tree->root;
		}
	}
	node->red = 0;

	retn:
	free (orig);
}

int rbforest_mk (struct rbforest *forest, unsigned n) {
	if (!(forest->trees = malloc (n * sizeof (struct rbtree)))) return 0;
	for (unsigned i = 0; i < n; ++i) rbtree_mk (forest->trees + i);
	forest->score = 0;
	forest->cnt = 0;
	forest->n = n;
	return 1;
}

void rbforest_fin (struct rbforest *forest) {
	for (unsigned i = 0; i < forest->n; ++i) rbtree_clear (forest->trees + i);
	free (forest->trees);
	forest->trees = NULL;
	forest->cnt = 0;
	forest->score = 0;
}

int rbforest_add (
	struct rbforest *forest, unsigned x1, unsigned x2, signed score
) {
	struct rbnode *node1, *node2;
	if (!(node1 = malloc (sizeof (struct rbnode)))) goto node1_err;
	if (!(node2 = malloc (sizeof (struct rbnode)))) goto node2_err;

	node1->peer = x2;
	node2->peer = x1;
	node1->cross = node2;
	node2->cross = node1;
	node1->score = node2->score = score;
	rbtree_add (forest->trees + x1, node1);
	rbtree_add (forest->trees + x2, node2);
	forest->score += score;
	--forest->cnt;

	return 1; node2_err:
	free (node1); node1_err:
	return 0;
}

void rbforest_rm (struct rbforest *forest, unsigned x) {
	struct rbtree *tree = forest->trees + x;
	for (
		struct rbnode *node = rbtree_first (tree);
		node; node = rbtree_next (tree, node)
	) {
		rbtree_rm (forest->trees + node->peer, node->cross);
		forest->score -= node->score;
		--forest->cnt;
	}
	rbtree_clear (tree);
}

void rbforest_pairs (struct rbforest const *forest, unsigned pairs[][2]) {
	struct rbtree *tree = forest->trees;
	unsigned idx = 0;
	for (unsigned i = 1; i < forest->n; ++i, ++tree) {
		for (
			struct rbnode *node = rbtree_first (tree);
			node && node->peer < i; node = rbtree_next (tree, node)
		) {
			pairs[idx][0] = i;
			pairs[idx++][1] = node->peer;
		}
	}
}

