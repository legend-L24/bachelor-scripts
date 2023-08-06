struct brin_t {
	struct crin_t *crin;
	float *rads;
};
typedef int narrow_f (crystal const *, unsigned, unsigned, signed *);
typedef int broad_f (crystal const *, narrow_f *, signed *);

extern struct bump_t *bump_mk (unsigned axes, unsigned n);
extern void bump_fin (struct bump_t *state);
extern int cryst_bmk (
	crystal *cr, float const rads[],
	struct peak_t const peaks[], struct atom_t const atoms[],
	float const abc[3], float const angles[3], unsigned nElem
);
extern void bryst_fin (crystal *cr);
extern struct brin_t *brin_read (FILE *file);
extern void brin_fin (struct brin_t *brin);
extern crystal *brin_eval (struct brin_t const *brin, rng_t *rng);

extern int bump_eval_orig (crystal const *cr, narrow_f *narrow, signed *ret);
extern int bump_eval_bc (crystal const *cr, narrow_f *narrow, signed *ret);
extern int bump_eval_bb (crystal const *cr, narrow_f *narrow, signed *ret);
extern int narrow_orig
	(crystal const *cr, unsigned i, unsigned j, signed *ret);
extern int narrow_cryst
	(crystal const *cr, unsigned i, unsigned j, signed *ret);
extern int cryst_eval_bb (crystal *cr, broad_f *broad, narrow_f *narrow);
extern int cryst_eval_bc (crystal *cr);

extern unsigned cryst_nball (crystal const *cr);
extern void cryst_set (crystal *cr, float const xyzs[][3]);
extern void cryst_write (crystal const *cr);

