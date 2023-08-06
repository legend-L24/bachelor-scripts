#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "rbtree.h"
#include "metric.h"
#include "utils.h"
#include "cryst.h"
#ifdef BENCHMARK
#include "bench_metric.h"
#include "bench_cryst.h"
#endif
#include "int_cryst.h"

struct cr_wptr {
	struct affine *affs;
	unsigned (*aids)[3];
};

static void *cp_alloc (void const *src, unsigned n) {
	void *buf;
	if (!(buf = malloc (n))) return NULL;
	memcpy (buf, src, n);
	return buf;
}

float torus_pos (float x) {
	return x > 0.0 ? x - (signed) x : x + (signed) (-x) + 1;
}

static float kx_dot (signed const k[3], float const x[3]) {
	float ret = 0.0;
	for (unsigned i = 0; i < 3; ++i) ret += k[i] * x[i];
	return 2.0 * PI * ret;
}

static float aff_dot (struct affine const *aff, float const vec[3]) {
	float ret = aff->move;
	struct linear const *lin = aff->lin;
	for (unsigned i = 0; i < aff->n; ++i) {
		ret += vec[lin[i].axis] * lin[i].zoom;
	}
	return ret;
}

static float orig_weight (
	signed const hkl[3], float const origs[][3], float nOrig
) {
	float ret[2] = { 0.0, 0.0 };
	for (unsigned i = 0; i < nOrig; ++i) {
		float kx = kx_dot (hkl, origs[i]);
		ret[0] += cos (kx);
		ret[1] += sin (kx);
	}
	return ret[0] * ret[0] + ret[1] * ret[1];
}

void balls_eval (
	struct crystal *cr, struct cr_atom *atom, unsigned nn[2]
) {
	struct cr_scat *scats = atom->scats;
	struct affine *affs = atom->affs;
	float (*origs)[3] = cr->origs;
	/*(cr->atoms)->xyz[0] = 0.0899484;
	(cr->atoms)->xyz[1] = 0.0899484;
	(cr->atoms)->xyz[2] = 0.156968;
        0.0102591       0.0102591       0.096851
	cr->atoms->xyz[0] = 0.741763;
        (cr->atoms+1)->xyz[0] = 0.80918;
        (cr->atoms+2)->xyz[0] = 0.264799;
        (cr->atoms+3)->xyz[1] = 0.551461;*/
	unsigned (*aids)[3] = atom->aids, i, j, k;
	nn[0] = atom->wyck[2] + 1;
	nn[1] = nn[0] * cr->nOrig;
	for (i = 0; i < atom->wyck[1]; ++i) {
		struct cr_ball *balls = scats[i].balls;
		for (j = 0; j < 3; ++j) {
			balls->xyz[j] = torus_pos (aff_dot (affs + aids[i][j], atom->xyz));
			if (cr->scales[2] == 0.0) continue;
			struct cr_ball *ball = balls;
			for (k = 1; k < cr->nOrig; ++k) {
				ball += nn[0];
				ball->xyz[j] = torus_pos (balls->xyz[j] + origs[k][j]);
			}
			if (atom->wyck[2]) {
				ball = balls;
				for (k = 0; k < cr->nOrig; ++k, ball += 2) {
					ball[1].xyz[j] = torus_pos (-ball[0].xyz[j]);
				}
			}
		}
	}
}
//this is to evaluate the scatter intensity of crystal surface
void surf_eval (struct crystal *cr, struct cr_atom *atom){
	struct cryst_surface *surf = cr->surf;
	unsigned nn[2], i, j;
	balls_eval (cr, atom, nn);
	for (i = 0; i < atom->wyck[1]; ++i) {
		struct cr_scat *scat = atom->scats + i;
		struct cr_ball *balls = scat->balls;
		for (j = 0; j < nn[1]; ++j) balls[j].flag |= 0x1;
		for (j = 0; j < cr->nSurf; ++j) if(cr->surf[j].flag){
			struct cryst_surface *surf_t = surf + j;
			float sasf = cr->scales[0] * atom->asfs[j],
				kx = kx_dot (surf_t->hkl, balls->xyz);
			signed *psc = surf_t->sc, *bsc = scat->scs[j];
			if (atom->wyck[2]) {
				psc[0] -= bsc[0];
				bsc[0] = 2 * sasf * cos (kx);
				psc[0] += bsc[0];
			} else {
				psc[0] -= bsc[0];
				psc[1] -= bsc[1];
				bsc[0] = sasf * cos (kx);
				bsc[1] = sasf * sin (kx);
				psc[0] += bsc[0];
				psc[1] += bsc[1];
			}
		}
	}
}

void atom_eval (struct crystal *cr, struct cr_atom *atom) {
	struct cr_peak *peaks = cr->hills->peaks;
	unsigned nn[2], i, j;
	balls_eval (cr, atom, nn);
	for (i = 0; i < atom->wyck[1]; ++i) {
		struct cr_scat *scat = atom->scats + i;
		struct cr_ball *balls = scat->balls;
		for (j = 0; j < nn[1]; ++j) balls[j].flag |= 0x1;
		for (j = 0; j < cr->nHill[1]; ++j) {
			struct cr_peak *peak = peaks + j;
			float sasf = cr->scales[0] * atom->asfs[j],
				kx = kx_dot (peak->peak.hkl, balls->xyz);
			signed *psc = peak->sc, *bsc = scat->scs[j];
			if (atom->wyck[2]) {
				psc[0] -= bsc[0];
				bsc[0] = 2 * sasf * cos (kx);
				psc[0] += bsc[0];
			} else {
				psc[0] -= bsc[0];
				psc[1] -= bsc[1];
				bsc[0] = sasf * cos (kx);
				bsc[1] = sasf * sin (kx);
				psc[0] += bsc[0];
				psc[1] += bsc[1];
			}
		}
	}
}

#ifdef BENCHMARK
void atom_eval_bb (struct crystal *cr, struct cr_atom *atom, int sc) {
	struct cr_peak *peaks = cr->hills->peaks;
	struct cr_scat *scats = atom->scats;
	unsigned nn[2], i, j, k;
	balls_eval (cr, atom, nn);
	for (i = 0; i < atom->wyck[1]; ++i) {
		struct cr_ball *balls = scats[i].balls;
		if (sc) for (j = 0; j < cr->nHill[1]; ++j) {
			struct cr_peak *peak = peaks + j;
			float asf = atom->asfs[j],
				kx = kx_dot (peak->peak.hkl, balls->xyz);
			if (atom->wyck[2]) peak->bsc[0] += 2 * asf * cos (kx);
			else {
				peak->bsc[0] += asf * cos (kx);
				peak->bsc[1] += asf * sin (kx);
			}
		} else for (j = 0; j < nn[1]; ++j) {
			struct cr_ball *ball = balls + j;
			for (k = 0; k < 3; ++k) {
				ball->bound[k][0] = torus_pos (ball->xyz[k] - atom->wid[k]);
				ball->bound[k][1] = torus_pos (ball->xyz[k] + atom->wid[k]);
			}
		}
	}
}
#endif

static struct crystal *cr_self_mk (
	struct hill_t const hills[], float const origs[][3],
	unsigned const wycks[][3], struct atom_t const atoms[], struct poly_t const polys[],
	float const abc[3], float const angles[3], float const ctl[3],
	rng_t *rng, unsigned nHill, unsigned nOrig, unsigned nAtom, unsigned nSpec, unsigned nSurf, unsigned nPoly
) {
	struct crystal *cr;
	unsigned i, j;
	if (!(cr = malloc (sizeof (struct crystal)))) goto cr_err;

	*cr = (struct crystal) {
		.nHill = { nHill, 0 }, .nAtom = { nAtom, 0 }, .nPoly = { nPoly, 0, 0, 0},
		.ctl = { ctl[1], ctl[2], 1.0 / (ctl[2] - ctl[1]) },
		.scores = { 0.0 }, .nOrig = nOrig, .rng = rng, .num = nSpec, .nSurf = nSurf// change something
	};
	for (i = 0; i < nHill; ++i) cr->nHill[1] += hills[i].n;
	for (cr->nPos[0] = i = 0; i < nAtom; ++i) for (j = 0; j < 3; ++j) {
		if (atoms[i].bound[j][0] != atoms[i].bound[j][1]) ++cr->nPos[0];
	}
        for (i = 0; i < nAtom; ++i) for (j = 0; j < nPoly; ++j) {
		if (atoms[i].elem == polys[j].center) {
                        unsigned const *wyck = wycks[polys[j].wyck];
			cr->nPoly[1] += 1;
                        cr->nPoly[2] += polys[j].n;
                        cr->nPoly[3] += polys[j].n * nOrig * wyck[1] * (wyck[2] + 1);
		}
	}
	for (i = 0; i < nAtom; ++i) {
		unsigned const *wyck = wycks[atoms[i].wyck];
		cr->nAtom[1] += nOrig * wyck[1] * (wyck[2] + 1);
	}
	cr->scales[1] = 0.99 * INT_MAX / ((cr->nAtom[1]+cr->nPoly[3]) * (cr->nAtom[1] + cr->nPoly[3] - 1));
	cr->scales[2] = ctl[0];

	if (!rbforest_mk (&cr->bump, cr->nAtom[1]+cr->nPoly[3])) goto bump_err;
	if (!(
		cr->poss = calloc (cr->nPos[0], sizeof (struct cr_pos))
	)) goto poss_err;
        if (!(
		cr->rots = calloc (3*cr->nPoly[1], sizeof (struct cr_rot))
	)) goto rots_err;
	if (!(
		cr->perm = calloc (cr->nPos[0]+3*cr->nPoly[1], sizeof (unsigned))
	)) goto perm_err;
	if (!(
		cr->origs = cp_alloc (origs, nOrig * 3 * sizeof (float))
	)) goto origs_err;
	if (!(cr->met = metric_mk (abc, angles))) goto met_err;

	for (i = 0; i < cr->nPos[0]+3*cr->nPoly[1]; ++i) cr->perm[i] = i;
	if (cr->nPos[0]) rng_shuf (rng, cr->perm, cr->nPos[0]+3*cr->nPoly[1]);
	cr->nPos[1] = 0;

	return cr;
	free (cr->met); met_err:
	free (cr->origs); origs_err:
	free (cr->perm); perm_err:
        free (cr->rots); rots_err:
	free (cr->poss); poss_err:
	rbforest_fin (&cr->bump); bump_err:
	free (cr); cr_err:
	return NULL;
}

static void cr_self_fin (struct crystal *cr) {
	free (cr->met);
	free (cr->origs);
	free (cr->perm);
	free (cr->poss);
	rbforest_fin (&cr->bump);
	free (cr);
}

static int cr_hills_mk (
	struct crystal *cr, struct peak_t const peaks[],
	struct hill_t const hills[]
) {
	if (!(
		cr->hills = calloc (cr->nHill[0], sizeof (struct cr_hill))
	)) goto hills_err;
	if (!(
		cr->hills->peaks = calloc (cr->nHill[1], sizeof (struct cr_peak))
	)) goto peaks_err;

	struct cr_peak *crPeaks = cr->hills->peaks;
	float sum = 0.0;
	unsigned i;

	for (i = 0; i < cr->nHill[0]; ++i) {
		if (i) cr->hills[i].peaks = cr->hills[i - 1].peaks + hills[i - 1].n;
		cr->hills[i].hill = hills[i];
		sum += hills[i].val0;
	}
	for (i = 0; i < cr->nHill[0]; ++i) cr->hills[i].hill.val0 /= sum;
	for (i = 0; i < cr->nHill[1]; ++i) {
		crPeaks[i].peak = peaks[i];
		crPeaks[i].peak.weight *= orig_weight
			(crPeaks[i].peak.hkl, cr->origs, cr->nOrig);
	}

	return 1;
	free (cr->hills->peaks); peaks_err:
	free (cr->hills); hills_err:
	return 0;
}

static void cr_hills_fin (struct crystal *cr) {
	free (cr->hills->peaks);
	free (cr->hills);
}
// lyt add in 2021/11/1
static int cr_surf_mk (
	struct crystal *cr, struct cryst_surface_t const surf[]
) {	unsigned i;
	if (!(
		cr->surf = calloc (cr->nSurf, sizeof(struct cryst_surface))
	)) goto surf_err;

	for (i = 0; i < cr->nSurf; ++i) {
		memcpy(cr->surf[i].hkl,surf[i].hkl,3*sizeof(signed int));
		cr->surf[i].mutil = surf[i].mutil * orig_weight
			(surf[i].hkl, cr->origs, cr->nOrig);
                // mutil=multiplicity * LP factor * origin_weight
		cr->surf[i].theta = surf[i].theta;
		if (!(
            		cr->surf[i].pos = calloc (1, sizeof (struct peak_shape))
        	)) goto pos_err;
	}
	return 1;
	free (cr->surf->pos); pos_err:
	free (cr->surf); surf_err:
	return 0;
}
// add the function about spectrum
// lyt add in 2021/10/15
static int cr_spec_mk (
	struct crystal *cr, struct spec_t const spec[]
) {
	float max = 0.0;
	unsigned i;
	if (!(
		cr->height = calloc (cr->num, sizeof(float))
	)) goto height_err;
	if (!(
		cr->cal_height = calloc (cr->num + 2*Val_spec, sizeof(float))
	)) goto cal_height_err;
	cr->cal_height += Val_spec;

	cr->origin = spec[0].theta;
        cr->inter = spec[1].theta-spec[0].theta;
	for (i = 0; i < cr->num; ++i) max = (spec[i].I > max) ? spec[i].I : max ;
	for (i = 0; i < cr->num; ++i) cr->height[i] = spec[i].I/max;
	return 1;
	free (cr->cal_height); cal_height_err:
	free (cr->height); height_err:
	return 0;
}

static void cr_surf_fin (struct crystal *cr) {
	free (cr->surf);
}

// add new function
static void cr_spec_fin (struct crystal *cr) {
	free (cr->height);
}

static int cr_atoms_mk (
	struct crystal *cr, struct affine const affs[], unsigned const aids[][3],
	unsigned const wycks[][3], float const dists[], float const asfs[], float const asfs_para[],
	struct atom_t const atoms[],struct poly_t const polys[],  unsigned nWyck, unsigned nElem
) {
	struct cr_scat *scats;
	unsigned i, j, k;

	for (j = i = 0; i < nWyck; ++i) j += wycks[i][0];
	if (!(
		cr->affs = cp_alloc (affs, j * sizeof (struct affine))
	)) goto affs_err;
	for (j = i = 0; i < nWyck; ++i) j += wycks[i][1];
	if (!(cr->aids = cp_alloc (aids, j * 3 * sizeof (unsigned)))) goto aids_err;

	if (!(
		cr->dists = cp_alloc (dists, nElem * nElem * sizeof (float))
	)) goto dists_err;
	if (!(cr->asfs = cp_alloc (
		asfs, nElem * cr->nSurf * sizeof (float)
	))) goto asfs_err; //lyt change nHill[1] to cr->nSurf

        if (!(cr->asfs_para = cp_alloc (
		asfs_para, nElem * 9 * sizeof (float)
	))) goto asfs_para_err; 
        cr->nElem = nElem;
	if (!(
		cr->atoms = calloc (cr->nAtom[0]+cr->nPoly[2], sizeof (struct cr_atom))
	)) goto atoms_err;
        if (!(
		cr->polys = calloc (cr->nPoly[1], sizeof (struct cr_poly))
	)) goto polys_err;
	for (j = i = 0; i < cr->nAtom[0]; ++i) j += wycks[atoms[i].wyck][1];
	for (i = 0; i < cr->nAtom[0]; ++i) for (k = 0; k < cr->nPoly[0]; ++k) {
		if (atoms[i].elem == polys[k].center) {
			j += polys[k].n * wycks[polys[k].wyck][1];
		}
	}
        if (!(
		scats = cr->atoms->scats = calloc (j, sizeof (struct cr_scat))
	)) goto scats_err;

	if (!(scats->scs = calloc (
		j * cr->nSurf, 2 * sizeof (signed)
	))) goto scs_err; //I change cr->nHill[1] to cr->nSurf I have not do everyone
	if (!(
		scats->balls = calloc (cr->nAtom[1]+cr->nPoly[3], sizeof (struct cr_ball))
	)) goto balls_err;

	return 1;
	free (cr->atoms->scats->balls); balls_err:
	free (cr->atoms->scats->scs); scs_err:
	free (cr->atoms->scats); scats_err:
        free (cr->polys); polys_err:
	free (cr->atoms); atoms_err:
        free (cr->asfs_para); asfs_para_err:
	free (cr->asfs); asfs_err:
	free (cr->dists); dists_err:
	free (cr->aids); aids_err:
	free (cr->affs); affs_err:
	return 0;
}

static void cr_atoms_fin (struct crystal *cr) {
	free (cr->atoms->scats->balls);
	free (cr->atoms->scats->scs);
	free (cr->atoms->scats);
	free (cr->atoms);
	free (cr->asfs);
	free (cr->dists);
	free (cr->aids);
	free (cr->affs);
}

static int cr_atoms_init (
	struct crystal *cr, unsigned const wycks[][3], struct atom_t const atoms[], 
        struct poly_t const polys[], struct extra_atom_t const extra_atoms[],
	unsigned nWyck, unsigned nElem
) {
	struct cr_atom *crAtoms = cr->atoms;
	struct cr_wptr *wps;
	unsigned i, j, k, m, n, q=0;

	if (!(wps = calloc (nWyck, sizeof (struct cr_wptr)))) return 0;
	wps[0] = (struct cr_wptr) { .affs = cr->affs, .aids = cr->aids };
	for (i = 1; i < nWyck; ++i) wps[i] = (struct cr_wptr) {
		.affs = wps[i - 1].affs + wycks[i - 1][0],
		.aids = wps[i - 1].aids + wycks[i - 1][1]
	};

	for (i = 0; i < cr->nAtom[0]; ++i) {
		struct cr_atom *atom = crAtoms + i;
		struct cr_scat *scat;
		unsigned elem = atoms[i].elem, wyck = atoms[i].wyck;
		if (i) {
			struct cr_atom *prev = atom - 1;
			scat = atom->scats = prev->scats + prev->wyck[1];
			scat->scs = prev->scats->scs + prev->wyck[1] * cr->nSurf;// // I change cr->nHill[1]
			scat->balls = prev->scats->balls + k;
		}
		atom->asfs = cr->asfs + elem * cr->nSurf; //lyt change nHill[1] to cr->nSurf
		memcpy (atom->wyck, wycks[wyck], 3 * sizeof (unsigned));
		atom->affs = wps[wyck].affs;
		atom->aids = wps[wyck].aids;
		k = cr->nOrig * (atom->wyck[2] + 1);
		for (j = 1; j < atom->wyck[1]; ++j) {
			scat = atom->scats + j;
			scat->scs = (scat - 1)->scs + cr->nSurf; // I change cr->nHill[1]
			scat->balls = (scat - 1)->balls + k;
		}

		scat = atom->scats;
		k *= atom->wyck[1];
		for (j = 0; j < k; ++j) scat->balls[j] = (struct cr_ball)
			{ .dists = cr->dists + elem * nElem, .elem = elem };
		scat->balls->n = k;
		scat->balls->flag = 0x2;
	}
        unsigned n_atom = 0;
	struct cr_atom *atom = crAtoms + cr->nAtom[0];
	struct cr_poly *poly = cr->polys;
	for (m = 0; m < cr->nPoly[0]; ++m) {
		unsigned elem = polys[m].center, wyck = polys[m].wyck, num = polys[m].n;
		for (n = 0; n < cr->nAtom[0]; ++n){
			if (elem == atoms[n].elem)	{
				poly->n = num;
				poly->orig = cr->atoms + n;
				poly->atom = atom;
				poly->ref_coor = calloc(num, 3 * sizeof (float));
				unsigned flag = 0;
				unsigned elem_0;
				for (i = 0; i < num; ++i) {
					memcpy (poly->ref_coor+i, (extra_atoms+i+n_atom)->xyz, 3 * sizeof (float));
					elem_0 = (extra_atoms+i+n_atom)->elem;
					struct cr_atom *prev = atom - 1;
					struct cr_scat *scat;
					scat = atom->scats = prev->scats + prev->wyck[1];
					scat->scs = prev->scats->scs + prev->wyck[1] * cr->nSurf;// // I change cr->nHill[1]
					scat->balls = prev->scats->balls + k;
					atom->asfs = cr->asfs + elem_0 * cr->nSurf; //lyt change nHill[1] to cr->nSurf
					memcpy (atom->wyck, wycks[wyck], 3 * sizeof (unsigned));
					atom->affs = wps[wyck].affs;
					atom->aids = wps[wyck].aids;
					k = cr->nOrig * (atom->wyck[2] + 1);
					for (j = 1; j < atom->wyck[1]; ++j) {
						scat = atom->scats + j;
						scat->scs = (scat - 1)->scs + cr->nSurf; // I change cr->nHill[1]
						scat->balls = (scat - 1)->balls + k;
					}
					scat = atom->scats;
					k *= atom->wyck[1];
					for (j = 0; j < k; ++j) scat->balls[j] = (struct cr_ball)
						{ .dists = cr->dists + elem_0 * nElem, .elem = elem_0 };
					scat->balls->n = k;
					scat->balls->flag = 0x2;
					++atom;
				}
				for (unsigned p = 0; p < 3; ++p) {
					if (polys[m].bound[p][0] == polys[m].bound[p][1]) {
						poly->angles[p] = polys[m].bound[p][0];
					} else{
						struct cr_rot *rot = cr->rots + q;
						*rot = (struct cr_rot) {
							.poly = poly, .axis = p, .flag = 1, .bound = {
								polys[m].bound[p][0], polys[m].bound[p][1],
								polys[m].bound[p][1] - polys[m].bound[p][0]
							}
						};
						rot->rot[0] = rot->rot[1] =
							rot->bound[0] + rot->bound[2] * rng_float (cr->rng);
						q++;
						flag=1; 
					}
				}
				if (!flag) {
					struct cr_atom *atom_0 = poly->atom;
					rot_mat_mk(poly->trans, cr->cell->para, cr->cell->para+3, poly->angles);
					for (unsigned r = 0; r < poly->n; ++r, ++atom_0) {
						mat_mult_add(poly->trans, poly->ref_coor[r], poly->orig->xyz, atom_0->xyz);
						atom_0->flag = 1;
					}
				}
				++poly;
			}
		}
		n_atom += polys[m].n;
	}
        poly = cr->polys;
	for (m = 0; m < cr->nPoly[1]; ++m) {
		unsigned  num = poly->n;
		struct cr_atom *atom = poly->orig;
		struct cr_scat *scats = atom->scats;
		for (n = 0; n < atom->wyck[1]; ++n){
			scats->balls[n].ind = m*atom->wyck[1]+n+1;
		}
		for (unsigned p = 0; p < num; ++p){
			struct cr_atom *atom = poly->atom + p;
			struct cr_scat *scats = atom->scats;
			for (n = 0; n < atom->wyck[1]; ++n){
				scats->balls[n].ind = m*atom->wyck[1]+n + 1;
			}
		}
		poly++;
	}
	free (wps);
	return 1;
}

static void cr_scale0_init (struct crystal *cr) {
	float sum = 0.0;
	for (unsigned i = 0; i < cr->nSurf; ++i) {  //lyt change nHill[1] to nSurf
		float s = 0.0;
		for (unsigned j = 0; j < cr->nAtom[0]+cr->nPoly[2]; ++j) {
			unsigned *wyck = cr->atoms[j].wyck;
			s += cr->atoms[j].asfs[i] * wyck[1] * (wyck[2] + 1);
		}
		if (sum < s) sum = s;
	}
        //cr->scales[0] = 1000;
	cr->scales[0] = 0.99 * INT_MAX / sum;
}

static void cr_poss_init (struct crystal *cr, struct atom_t const atoms[]) {
	for (unsigned i = 0, j = 0; j < cr->nAtom[0]; ++j) {
		struct atom_t const *atom = atoms + j;
		struct cr_atom *crAtom = cr->atoms + j;
		int flag = 0;
		for (unsigned k = 0; k < 3; ++k) {
			if (atom->bound[k][0] == atom->bound[k][1]) {
				crAtom->xyz[k] = atom->bound[k][0];
			} else {
				struct cr_pos *pos = cr->poss + i;
				*pos = (struct cr_pos) {
					.atom = crAtom, .axis = k, .flag = 1, .bound = {
						atom->bound[k][0], atom->bound[k][1],
						atom->bound[k][1] - atom->bound[k][0]
					}
				};
				pos->pos[0] = pos->pos[1] =
					pos->bound[0] + pos->bound[2] * rng_float (cr->rng);
				flag = 1;
				++i;
			}
		}
                if (!flag) surf_eval(cr, crAtom); //atom_eval (cr, crAtom);
	}
}
static void cr_vari_fin (struct crystal *cr) {
	free (cr->peak_para->var);
	free (cr->peak_para->para);
	free (cr->cell->para);
	free (cr->cell->lin);
	free (cr->cell->var);
	free (cr->cell);
	free (cr->peak_para);
}
// U need to add boundary for varibles
static int cr_vari_init (struct crystal *cr, float const abc[3], float const angles[3], unsigned sg_number){
	unsigned i;
	cr->peak_para = calloc(1, sizeof(struct cr_peak_para));
	cr->cell = calloc(1, sizeof(struct cr_cell));
	cr->peak_para->n = 5;//U, V, W, X, Y are five parameters of cw2_fun
	cr->peak_para->var = calloc(cr->peak_para->n, sizeof(struct cr_var));
	cr->peak_para->para = calloc(cr->peak_para->n, sizeof(float));
	cr->cell->para = calloc(6, sizeof(float));
	memcpy (cr->cell->para,abc,sizeof (float)*3);
	memcpy (cr->cell->para+3,angles,sizeof(float)*3);
	// I will change 1 to 0 in the near future
	if (sg_number>=231) 
		return 0;
	else if (sg_number>=195)
		{	cr->cell->n = 1;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[1].axis = 1;
			cr->cell->lin[2].axis = 1;
			cr->cell->lin[3].zoom = 90;
			cr->cell->lin[4].zoom = 90;
			cr->cell->lin[5].zoom = 90;
			cr->cell->var->bound[0] = 15; cr->cell->var->bound[1] = 60; cr->cell->var->bound[2] = 45;
		}
	else if (sg_number>=168)
		{	cr->cell->n = 2;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = abc[2];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[1].axis = 1;
			cr->cell->lin[8].axis = 1;
			cr->cell->lin[3].zoom = 90;
			cr->cell->lin[4].zoom = 90;
			cr->cell->lin[5].zoom = 120;
		}
	else if (sg_number>=143)
		{	cr->cell->n = 2;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = angles[0];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[1].axis = 1;
			cr->cell->lin[2].axis = 1;
			cr->cell->lin[9].axis = 1;
			cr->cell->lin[10].axis = 1;
			cr->cell->lin[11].axis = 1;
		}
	else if (sg_number>=75)
		{	cr->cell->n = 2;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = abc[2];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[1].axis = 1;
			cr->cell->lin[8].axis = 1;
			cr->cell->lin[3].zoom = 90;
			cr->cell->lin[4].zoom = 90;
			cr->cell->lin[5].zoom = 90;
		}
	else if (sg_number>=16)
		{	cr->cell->n = 3;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = abc[1];
			(cr->cell->var+2)->val[0] = (cr->cell->var+2)->val[1] = abc[2];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[7].axis = 1;
			cr->cell->lin[14].axis = 1;
			cr->cell->lin[3].zoom = 90;
			cr->cell->lin[4].zoom = 90;
			cr->cell->lin[5].zoom = 90;
		}
	else if (sg_number>=3)
		{	cr->cell->n = 4;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = abc[1];
			(cr->cell->var+2)->val[0] = (cr->cell->var+2)->val[1] = abc[2];
			(cr->cell->var+3)->val[0] = (cr->cell->var+3)->val[1] = angles[1];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[7].axis = 1;
			cr->cell->lin[14].axis = 1;
			cr->cell->lin[3].zoom = 90;
			cr->cell->lin[22].axis = 1;
			cr->cell->lin[5].zoom = 90;
		}
	else if (sg_number>=1)
		{	cr->cell->n = 6;
			cr->cell->var = calloc(cr->cell->n, sizeof(struct cr_var));
			cr->cell->lin = calloc(cr->cell->n*6, sizeof(struct linear));
			cr->cell->var->val[0] = cr->cell->var->val[1] = abc[0];
			(cr->cell->var+1)->val[0] = (cr->cell->var+1)->val[1] = abc[1];
			(cr->cell->var+2)->val[0] = (cr->cell->var+2)->val[1] = abc[2];
			(cr->cell->var+3)->val[0] = (cr->cell->var+3)->val[1] = angles[0];
			(cr->cell->var+4)->val[0] = (cr->cell->var+4)->val[1] = angles[1];
			(cr->cell->var+5)->val[0] = (cr->cell->var+5)->val[1] = angles[2];
			cr->cell->lin[0].axis = 1;
			cr->cell->lin[7].axis = 1;
			cr->cell->lin[14].axis = 1;
			cr->cell->lin[21].axis = 1;
			cr->cell->lin[28].axis = 1;
			cr->cell->lin[35].axis = 1;
		}
	else goto assign_err;

	for (i =0; i < cr->cell->n; ++i){
		cr->cell->var[i].flag = 1;
	}
	cr->nCell[0] = cr->cell->n; 
	cr->nCell[1] = 0;
	if (!(
		cr->vari_perm = calloc (cr->nCell[0], sizeof (unsigned))
	)) goto vari_perm_err;
	for (i = 0; i < cr->nCell[0]; ++i) cr->vari_perm[i] = i;
	if (cr->nCell[0]) rng_shuf (cr->rng, cr->vari_perm, cr->nCell[0]);
	return 1;

	free(cr->vari_perm); vari_perm_err:
	free(cr->cell); assign_err:
	return 0;
}

struct crystal *cryst_mk (
	struct peak_t const peaks[], struct hill_t const hills[], struct spec_t const spec[],
	float const origs[][3], struct affine const affs[],
	unsigned const aids[][3], unsigned const wycks[][3],
	float const dists[], float const asfs[], float const asfs_para[], struct atom_t const atoms[], 
        struct poly_t const polys[], struct extra_atom_t const extra_atoms[], struct cryst_surface_t const surf[],
	float const abc[3], float const angles[3], float const ctl[3],
	rng_t *rng, unsigned nHill, unsigned nOrig,
	unsigned nWyck, unsigned nElem, unsigned nAtom, unsigned nSpec, unsigned nSurf, unsigned sg_number, unsigned nPoly
) {
	struct crystal *cr;

	if (!(cr = cr_self_mk (
		hills, origs, wycks, atoms, polys, abc, angles,
		ctl, rng, nHill, nOrig, nAtom, nSpec, nSurf, nPoly
	))) goto cr_err;
	if (!cr_hills_mk (cr, peaks, hills)) goto hills_err;
	if (!cr_surf_mk (cr, surf)) goto surf_err; //change
        if (!cr_spec_mk (cr, spec)) goto spec_err;//add by lyt
	if (!cr_atoms_mk (
		cr, affs, aids, wycks, dists, asfs, asfs_para, atoms, polys, nWyck, nElem
	)) goto atoms_err;
	if (!cr_vari_init (cr, abc, angles, sg_number)) goto vari_init_err;
	if (!cr_atoms_init (cr, wycks, atoms, polys, extra_atoms, nWyck, nElem)) goto init_err;
	cr_scale0_init (cr);
	cr_poss_init (cr, atoms);
        return cr; init_err:
        cr_vari_fin(cr); vari_init_err:
	cr_atoms_fin (cr); atoms_err:
        // add cr_spec_fin
	cr_spec_fin (cr); spec_err:
	cr_surf_fin (cr); surf_err:
	cr_hills_fin (cr); hills_err:
	cr_self_fin (cr); cr_err:
	return NULL;
}

void cryst_fin (struct crystal *cr) {
	cr_atoms_fin (cr);
	cr_hills_fin (cr);
	cr_self_fin (cr);
}

#ifdef BENCHMARK
int cryst_bmk (
	struct crystal *cr, float const rads[],
	struct peak_t const peaks[], struct atom_t const atoms[],
	float const abc[3], float const angles[3], unsigned nElem
) {
	float ratios[3];
	unsigned axes = 0x0, i, j;

	ratios_get (abc, angles, ratios);
	for (j = 0, i = 1; i < 3; ++i) if (ratios[i] < ratios[j]) j = i;
	axes = 0x1 << j;
	if (!(cr->bbump = bump_mk (axes, cr->nAtom[1]))) goto bbump_err;
	if (!(cr->wids = calloc (nElem, 3 * sizeof (float)))) goto wids_err;

	for (i = 0; i < nElem; ++i) {
		for (j = 0; j < 3; ++j) if (
			(cr->wids[i][j] = cr->ctl[1] * ratios[j] * rads[i]) >= 0.4999
		) cr->wids[i][j] = 0.4999;
	}
	for (i = 0; i < cr->nAtom[0]; ++i) {
		cr->atoms[i].wid = cr->wids[atoms[i].elem];
	}

	return 1;
	free (cr->wids); wids_err:
	bump_fin (cr->bbump); bbump_err:
	return 0;
}

void bryst_fin (struct crystal *cr) {
	free (cr->wids);
	bump_fin (cr->bbump);
	cryst_fin (cr);
}
#endif

