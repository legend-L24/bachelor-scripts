#include <stdio.h>
#include <stdlib.h>
#include "rbtree.h"

#define LEN 10

static int arr_read (unsigned arr[], unsigned n) {
	for (unsigned i = 0; i < n; ++i) if (scanf ("%u", arr + i) != 1) return 0;
	return 1;
}

static void arr_write (unsigned const arr[], unsigned n) {
	if (!n) return;
	for (unsigned i = 0;;) {
		printf ("%u", arr[i]);
		if (++i < n) {
			printf (" ");
		} else {
			printf ("\n");
			break;
		}
	}
}

static int arr_chk (unsigned const arr[], unsigned n) {
	for (unsigned i = 1; i < n; ++i) if (arr[i - 1] >= arr[i]) return 0;
	return 1;
}

static int rbtree_dump (struct rbtree const *tree, unsigned arr[]) {
	unsigned i = 0;
	for (
		struct rbnode *node = rbtree_first (tree); node;
		node = rbtree_next (tree, node)
	) arr[i++] = node->peer;
	return i;
}

static int test_add (struct rbtree *tree, unsigned arr[3][LEN]) {
	for (unsigned i = 0; i < LEN; ++i) {
		struct rbnode *node;
		if (!(node = malloc (sizeof (struct rbnode)))) return 0;
		node->peer = arr[0][i];
		rbtree_add (tree, node);

		unsigned n = rbtree_dump (tree, arr[2]);
		if (!(n == i + 1 && arr_chk (arr[2], n))) {
			printf ("add:\n");
			arr_write (arr[0], i + 1);
			arr_write (arr[2], n);
			return 0;
		}
	}
	return 1;
}

static int test_rm (struct rbtree *tree, unsigned arr[3][LEN]) {
	for (unsigned i = 0; i < LEN; ++i) {
		struct rbnode *node;
		if (!(node = rbtree_get (tree, arr[1][i]))) return 0;
		rbtree_rm (tree, node);

		unsigned n = rbtree_dump (tree, arr[2]);
		if (!(n + i + 1 == LEN && arr_chk (arr[2], n))) {
			printf ("rm:\n");
			arr_write (arr[0], LEN);
			arr_write (arr[1], i);
			arr_write (arr[2], n);
			return 0;
		}
	}
	return 1;
}

int main (void) {
	struct rbtree tree;
	unsigned arr[3][LEN];

	rbtree_mk (&tree);
	return !(
		arr_read (arr[0], LEN) &&
		arr_read (arr[1], LEN) &&
		test_add (&tree, arr) &&
		test_rm (&tree, arr)
	);
}

