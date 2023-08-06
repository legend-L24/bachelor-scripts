#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
#include "utils.h"
#include "cryst.h"
#ifdef BENCHMARK
#include "bench_cryst.h"
#endif

static int crin_chk (struct crin_t const *crin, unsigned const n[3]) {
	unsigned i, j;
	int ret = 1;

	for (i = 0; i < 3; ++i) if (crin->abc[i] <= 0.0)
		{ fprintf (stderr, "E: invalid abc[%u]\n", i); ret = 0; }
	for (i = 0; i < 3; ++i) if (!(
		crin->angles[i] > 0.0 && crin->angles[i] < 180.0
	)) { fprintf (stderr, "E: invalid angles[%u]\n", i); ret = 0; }
	if (!(
		crin->ctl[0] >= 0.0 && crin->ctl[0] < 1.0 &&
		0.0 <= crin->ctl[1] && crin->ctl[1] < crin->ctl[2]
	)) { fprintf (stderr, "E: invalid `ctl'\n"); ret = 0; }
	for (i = 0; i < 3; ++i) if (crin->origs[0][i] != 0.0)
		{ fprintf (stderr, "E: invalid origs[0]\n"); ret = 0; }

	for (i = 0; i < crin->nHill; ++i) if (!(
		crin->hills[i].n && crin->hills[i].val0 >= 0.0
	)) { fprintf (stderr, "E: invalid hills[%u]\n", i); ret = 0; }
	for (i = 0; i < n[0]; ++i) if (crin->peaks[i].weight <= 0.0)
		{ fprintf (stderr, "E: invalid peaks[%u]\n", i); ret = 0; }
	for (i = 0; i < crin->nElem; ++i) {
		for (j = i; j < crin->nElem; ++j) {
			if (crin->dists[i * crin->nElem + j] <= 0.0) {
				fprintf (stderr, "E: invalid dists[%u][%u]\n", i, j); ret = 0;
			}
			if (
				crin->dists[i * crin->nElem + j] !=
					crin->dists[j * crin->nElem + i]
			) {
				fprintf (
					stderr, "E: dists[%u][%u] != dists[%u][%u]\n", i, j, j, i
				); ret = 0;
			}
		}
	}
	// lyt change the check program,n[0] to crin->nSurf
	for (i = 0; i < crin->nElem; ++i) for (j = 0; j < crin->nSurf; ++j) if (
		crin->asfs[i * crin->nSurf + j] <= 0.0
	) { fprintf (stderr, "E: invalid asfs[%u][%u]\n", i, j); ret = 0; }

	for (i = 0; i < crin->nAtom; ++i) {
		if (!(
			crin->atoms[i].elem < crin->nElem &&
			crin->atoms[i].wyck < crin->nWyck
		)) goto atoms_err;
		for (j = 0; j < 3; ++j) if (
			crin->atoms[i].bound[j][0] > crin->atoms[i].bound[j][1]
		) goto atoms_err;
		continue; atoms_err:
		fprintf (stderr, "E: invalid atoms[%u]\n", i); ret = 0;
	}
	for (i = 0; i < crin->nWyck; ++i) if (!(
		crin->wycks[i][0] && crin->wycks[i][1] && crin->wycks[i][2] < 2
	)) { fprintf (stderr, "E: invalid wycks[%u]\n", i); ret = 0; }
	for (i = 0; i < n[1]; ++i) {
		struct affine const *aff = crin->affs + i;
		unsigned mark[3] = { 0 };
		if (aff->n > 2) goto affs_err;
		for (j = 0; j < aff->n; ++j) {
			if (aff->lin[j].axis > 2 || mark[aff->lin[j].axis]) goto affs_err;
			mark[aff->lin[j].axis] = 1;
		}
		continue; affs_err:
		fprintf (stderr, "E: invalid affs[%u]\n", i); ret = 0;
	}
	for (j = i = 0; i < crin->nWyck; ++i) {
		for (unsigned end = j + crin->wycks[i][1]; j < end; ++j) {
			for (unsigned k = 0; k < 3; ++k) if (
				crin->aids[j][k] >= crin->wycks[i][0]
			) goto aids_err;
		}
		continue; aids_err:
		fprintf (stderr, "E: invalid aids[%u]\n", j); ret = 0;
	}

	return ret;
}

static int scan_orig (char const **ptr, float orig[3]) {
	return scan_farr (ptr, orig, 3);
}

static int scan_hill (char const **ptr, struct hill_t *hill) {
	return scan_unsigned (ptr, &hill->n) &&
		scan_space (ptr) && scan_float (ptr, &hill->val0);
}

static int scan_peak (char const **ptr, struct peak_t *peak) {
	return scan_iarr (ptr, peak->hkl, 3) &&
		scan_space (ptr) && scan_float (ptr, &peak->weight);
}
//lyt add 2021 11 1
static int scan_surf (char const **ptr, struct cryst_surface_t *surf) {
	return  scan_iarr (ptr, surf->hkl, 3)&&
		scan_space (ptr) && scan_float (ptr, &surf->mutil) && scan_float (ptr, &surf->theta);
}
//lyt add 2021 10 9
static int scan_spec (char const **ptr, struct spec_t *spec) {
	return scan_float (ptr, &spec->theta) &&
		scan_space (ptr) && scan_float (ptr, &spec->I);
}
static int scan_poly (char const **ptr, struct poly_t *polys) {
	return scan_unsigned (ptr, &polys->center) &&
		scan_space (ptr) && scan_unsigned (ptr, &polys->wyck) &&
		scan_space (ptr) && scan_unsigned (ptr, &polys->n) &&
		scan_space (ptr) && scan_farr (ptr, (float *) &polys->bound, 6);
}
static int scan_extra_atom (char const **ptr, struct extra_atom_t *extra_atoms) {
	return scan_unsigned (ptr, &extra_atoms->elem) &&
		scan_space (ptr) && scan_farr (ptr, (float *) &extra_atoms->xyz, 3);
}

static int scan_atom (char const **ptr, struct atom_t *atom) {
	return scan_unsigned (ptr, &atom->elem) &&
		scan_space (ptr) && scan_unsigned (ptr, &atom->wyck) &&
		scan_space (ptr) && scan_farr (ptr, (float *) &atom->bound, 6);
}

static int scan_u3 (char const **ptr, unsigned us[3]) {
	return scan_uarr (ptr, us, 3);
}

static int scan_aff (char const **ptr, struct affine *aff) {
	return scan_unsigned (ptr, &aff->n) &&
		scan_space (ptr) && scan_float (ptr, &aff->move) &&
		scan_space (ptr) && scan_unsigned (ptr, &aff->lin[0].axis) &&
		scan_space (ptr) && scan_float (ptr, &aff->lin[0].zoom) &&
		scan_space (ptr) && scan_unsigned (ptr, &aff->lin[1].axis) &&
		scan_space (ptr) && scan_float (ptr, &aff->lin[1].zoom);
}

static int crin_head_read (struct crin_t *crin, struct cursor *curs) {
	if (!cursor_step (curs)) return
		cursor_err (curs, "E: premature EOF"), 0;
	else if (scan_space (CUR_PTR (curs))) return
		cursor_err (curs, "E: leading whitespace"), 0;
	else if (!scan_farr (CUR_PTR (curs), crin->abc, 3)) return
		cursor_err (curs, "E: malformed `abc'"), 0;
	else if (!(scan_space (CUR_PTR (curs)) && scan_farr (
		CUR_PTR (curs), crin->angles, 3
	))) return cursor_err (curs, "E: malformed `angles'"), 0;
	else if (!(scan_space (CUR_PTR (curs)) && scan_farr (
		CUR_PTR (curs), crin->ctl, 3
	))) return cursor_err (curs, "E: malformed `ctl'"), 0;

	//lyt add in 2022/1/17
	else if (!(scan_space (CUR_PTR (curs)) && scan_unsigned (
		CUR_PTR (curs), &crin->sg_number
	))) return cursor_err (curs, "E: malformed `space group number'"), 0;

	else if (*curs->ptr) return
		cursor_err (curs, "E: trailing character"), 0;
	else return 1;
}

#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
#define SCAN_VAR(FIELD, SCAN, CNT, SIZE) \
	if (!(scan_block ( \
		&curs, SIZE, 0, &CNT, (void **) &crin->FIELD, (scan_f *) SCAN \
	) && CNT)) ERR_GOTO ("E: failed to parse `" #FIELD "'", FIELD##_err)
#define SCAN_FIX(FIELD, SCAN, CNT, SIZE) \
	if (!scan_block ( \
		&curs, SIZE, CNT, NULL, (void **) &crin->FIELD, (scan_f *) SCAN \
	)) ERR_GOTO ("E: failed to parse `" #FIELD "'", FIELD##_err)
struct crin_t *crin_read (FILE *file) {
	struct crin_t *crin;
	struct cursor curs;
	char const *desc = NULL;
	float nElem;
	unsigned n[4], i;

	cursor_mk (&curs, file);
	if (!(crin = malloc (sizeof (struct crin_t)))) {
		fprintf (stderr, "E: failed malloc() in %s()\n", __func__);
		goto crin_err;
	}
	if (!crin_head_read (crin, &curs)) ERR_GOTO
		("E: failed to parse header", origs_err)
	SCAN_VAR (origs, &scan_orig, crin->nOrig, 3 * sizeof (float))

	SCAN_VAR (hills, &scan_hill, crin->nHill, sizeof (struct hill_t))
	for (n[0] = i = 0; i < crin->nHill; ++i) n[0] += crin->hills[i].n;
	SCAN_FIX (peaks, &scan_peak, n[0], sizeof (struct peak_t))
        SCAN_VAR (surf, &scan_surf, crin->nSurf, sizeof (struct cryst_surface_t)) //change
	SCAN_VAR (spec, &scan_spec, crin->nSpec, sizeof (struct spec_t)) //change
	SCAN_VAR (dists, &scan_float, crin->nElem, sizeof (float))
	nElem = round (sqrt (crin->nElem));
	if (nElem * nElem == crin->nElem) crin->nElem = nElem;
	else ERR_GOTO ("E: invalid length of `dists'", retn_err)
	SCAN_FIX (asfs, &scan_float, crin->nElem * crin->nSurf, sizeof (float)) //lyt change n[0] to crin->nSurf
	SCAN_FIX (asfs_para, &scan_float, crin->nElem * 9, sizeof (float)) 
	SCAN_VAR (wycks, &scan_u3, crin->nWyck, 3 * sizeof (unsigned))
	for (n[1] = n[2] = i = 0; i < crin->nWyck; ++i) {
		n[1] += crin->wycks[i][0];
		n[2] += crin->wycks[i][1];
	}
        SCAN_FIX (affs, &scan_aff, n[1], sizeof (struct affine))
	SCAN_FIX (aids, &scan_u3, n[2], 3 * sizeof (unsigned))
	SCAN_VAR (atoms, &scan_atom, crin->nAtom, sizeof (struct atom_t))
        //SCAN_VAR (polys, &scan_poly, crin->nPoly, sizeof (struct poly_t))
	scan_block ( &curs, sizeof (struct poly_t), 0, &crin->nPoly, (void **) &crin->polys, (scan_f *) &scan_poly );
        for (n[3] = i = 0; i < crin->nPoly; ++i) {
		n[3] += crin->polys[i].n;
	}
        scan_block ( &curs, sizeof (struct extra_atom_t), n[3], NULL, (void **) &crin->extra_atoms, (scan_f *) &scan_extra_atom );
	//SCAN_FIX (extra_atoms, &scan_extra_atom, n[3], sizeof (struct extra_atom_t))
	if (cursor_step (&curs)) ERR_GOTO ("E: trailing line", retn_err)
	if (!crin_chk (crin, n)) goto retn_err;

	return crin; retn_err:
        //free (crin->extra_atoms); extra_atoms_err:
	//free (crin->polys); polys_err:
	free (crin->atoms); atoms_err:
	free (crin->aids); aids_err:
	free (crin->affs); affs_err:
	free (crin->wycks); wycks_err:
        free (crin->asfs_para); asfs_para_err:
	free (crin->asfs); asfs_err:
	free (crin->dists); dists_err:
        free (crin->spec); spec_err: // change
	free (crin->surf); surf_err: //change
	free (crin->peaks); peaks_err:
	free (crin->hills); hills_err:
	free (crin->origs); origs_err:
	free (crin); crin_err:
	if (desc) cursor_err (&curs, desc);
	return NULL;
}

void crin_fin (struct crin_t *crin) {
	free (crin->atoms);
	free (crin->aids);
	free (crin->affs);
	free (crin->wycks);
	free (crin->asfs);
	free (crin->dists);
	free (crin->peaks);
	free (crin->hills);
	free (crin->origs);
	free (crin);
}

crystal *crin_eval (struct crin_t const *crin, rng_t *rng) {
	return cryst_mk (
		crin->peaks, crin->hills, crin->spec, crin->origs,
		crin->affs, crin->aids, crin->wycks,
		crin->dists, crin->asfs, crin->asfs_para, crin->atoms, crin->polys, crin->extra_atoms, crin->surf,
		crin->abc, crin->angles, crin->ctl, rng,
		crin->nHill, crin->nOrig, crin->nWyck, crin->nElem, crin->nAtom, crin->nSpec, crin->nSurf, crin->sg_number, crin->nPoly
	);
}

crystal *cryst_read (rng_t *rng, FILE *file) {
	struct crin_t *crin = crin_read (file);
	crystal *cr = crin ? crin_eval (crin, rng) : NULL;
	if (crin) crin_fin (crin);
        return cr;
}

#ifdef BENCHMARK
struct brin_t *brin_read (FILE *file) {
	struct brin_t *brin;
	struct cursor curs;
	unsigned nElem;

	cursor_mk (&curs, file);
	if (!(brin = malloc (sizeof (struct brin_t)))) goto brin_err;
	if (!scan_block (
		&curs, sizeof (float), 0, &nElem,
		(void **) &brin->rads, (scan_f *) &scan_float
	)) goto rads_err;
	if (!(brin->crin = crin_read (file))) goto crin_err;
	if (brin->crin->nElem != nElem) goto retn_err;

	return brin; retn_err:
	crin_fin (brin->crin); crin_err:
	free (brin->rads); rads_err:
	free (brin); brin_err:
	return NULL;
}

void brin_fin (struct brin_t *brin) {
	crin_fin (brin->crin);
	free (brin->rads);
	free (brin);
}

crystal *brin_eval (struct brin_t const *brin, rng_t *rng) {
	crystal *cr;
	struct crin_t const *crin = brin->crin;
	if (!(cr = crin_eval (crin, rng))) goto cr_err;
	if (!cryst_bmk (
		cr, brin->rads, crin->peaks, crin->atoms,
		crin->abc, crin->angles, crin->nElem
	)) goto retn_err;
	return cr; retn_err:
	cryst_fin (cr); cr_err:
	return NULL;
}
#endif

