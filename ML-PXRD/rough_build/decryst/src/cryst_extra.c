#include <math.h>
#ifdef BENCHMARK
#include <stdlib.h>
#include <string.h>
#endif

#include "rng.h"
#include "rbtree.h"
#include "metric.h"
#include "utils.h"
#include "cryst.h"
#ifdef BENCHMARK
#include "bench_cryst.h"
#endif
#include "int_cryst.h"

#include <string.h>
#include <stdlib.h> // this is need by peak_fun



// I create some constants for profile
#define Pi 3.1415926
#define Ln2 0.6931

#ifdef BENCHMARK

struct sweep {
	float pos;
	unsigned idx;
	int end;
};

struct cpair {
	unsigned idx[2];
};

struct bump_t {
	struct rbtree active;
	struct sweep *poss;
	struct cpair *bumps[3];
	unsigned cnt[3], n;
};

struct bump_t *bump_mk (unsigned axes, unsigned n) {
	struct bump_t *state;
	unsigned i, j;

	if (!(state = malloc (sizeof (struct bump_t)))) goto state_err;
	for (i = 0; i < 3; ++i) state->bumps[i] = NULL;
	for (i = 0; i < 3; ++i) if (axes >> i & 0x1 && !(
		state->bumps[i] = calloc (n * n, sizeof (struct cpair))
	)) goto poss_err;
	if (!(state->poss = malloc (2 * n * sizeof (struct sweep)))) goto poss_err;

	for (j = i = 0; i < n; ++i) {
		state->poss[j++] = (struct sweep) { .idx = i, .end = 0 };
		state->poss[j++] = (struct sweep) { .idx = i, .end = 1 };
	}
	rbtree_mk (&state->active);
	state->n = n;
	return state;

	free (state->poss);
	poss_err:
	for (i = 0; i < 3; ++i) if (state->bumps[i]) free (state->bumps[i]); 
	free (state);
	state_err:
	return NULL;
}

void bump_fin (struct bump_t *state) {
	unsigned i;
	free (state->poss);
	for (i = 0; i < 3; ++i) if (state->bumps[i]) free (state->bumps[i]);
	free (state);
}

int bump_eval_orig (struct crystal const *cr, narrow_f *narrow, signed *ret_) {
	struct bump_t const *state = cr->bbump;
	signed ret = 0, score;
	for (unsigned i = 0; i < state->n; ++i) {
		for (unsigned j = i + 1; j < state->n; ++j) {
			if ((*narrow) (cr, i, j, &score)) ret += score;
		}
	}
	*ret_ = ret;
	return 1;
}

int bump_eval_bc (struct crystal const *cr, narrow_f *narrow, signed *ret_) {
	struct bump_t const *state = cr->bbump;
	struct cr_ball const *balls = cr->atoms->scats->balls;
	signed ret = 0, score;
	for (unsigned i = 0; i < state->n; ++i) if (balls[i].flag & 0x2) {
		for (unsigned j = 0; j < state->n; ++j) {
			if (i == j || (balls[j].flag & 0x2 && i > j)) continue;
			if ((*narrow) (cr, i, j, &score)) ret += score;
		}
	}
	*ret_ = ret;
	return 1;
}

static int sweep_cmp (void const *p1, void const *p2) {
	struct sweep const *s1 = p1, *s2 = p2;
	float const x = s1->pos - s2->pos;
	return x != 0.0 ? x > 0.0 : s2->end - s1->end;
}

static int cpair_cmp (void const *p1, void const *p2) {
	struct cpair const *c1 = p1, *c2 = p2;
	return c1->idx[0] != c2->idx[0] ?
		c1->idx[0] - c2->idx[0] : c1->idx[1] - c2->idx[1];
}

static int bump_axis_bb (
	struct bump_t *state, struct cr_ball const balls[], unsigned axis
) {
	struct rbtree *active = &state->active;
	struct rbnode *node, *tmp;
	struct cpair *bump = state->bumps[axis];
	unsigned nn = 2 * state->n, cnt, i, j;

	for (j = i = 0; i < state->n; ++i) {
		state->poss[j++] = (struct sweep)
			{ .idx = i, .end = 0, .pos = balls[i].bound[axis][0] };
		state->poss[j++] = (struct sweep)
			{ .idx = i, .end = 1, .pos = balls[i].bound[axis][1] };
	}
	qsort (state->poss, nn, sizeof (struct sweep), &sweep_cmp);
	cnt = 0;

	for (i = 0; i < nn; ++i) {
		j = state->poss[i].idx;
		if (state->poss[i].end) {
			if ((tmp = rbtree_get (active, j))) rbtree_rm (active, tmp);
		} else {
			if (!(tmp = malloc (sizeof (struct rbnode)))) return 0;
			for (
				node = rbtree_first (active); node;
				node = rbtree_next (active, node)
			) if ((balls[j].flag | balls[node->peer].flag) & 0x2) {
				bump[cnt++] = j < node->peer ?
					(struct cpair) {{ j, node->peer }} :
					(struct cpair) {{ node->peer, j }};
			}
			tmp->peer = j;
			rbtree_add (active, tmp);
		}
	}

	// Here duplicates might be introduced.
	for (i = 0; active->root != &active->nil; ++i) {
		j = state->poss[i].idx;
		if (state->poss[i].end) {
			if ((tmp = rbtree_get (active, j))) rbtree_rm (active, tmp);
		} else if (!rbtree_get (active, j)) for (
			node = rbtree_first (active); node;
			node = rbtree_next (active, node)
		) if ((balls[j].flag | balls[node->peer].flag) & 0x2) {
			bump[cnt++] = j < node->peer ?
				(struct cpair) {{ j, node->peer }} :
				(struct cpair) {{ node->peer, j }};
		}
	}

	rbtree_clear (active);
	qsort (bump, cnt, sizeof (struct cpair), &cpair_cmp);
	state->cnt[axis] = cnt;
	return 1;
}

int bump_eval_bb (struct crystal const *cr, narrow_f *narrow, signed *ret_) {
	struct bump_t *state = cr->bbump;
	struct cr_ball const *balls = cr->atoms->scats->balls;
	struct cpair *cur[3], *end[3], cp = {{ 0, 0 }};
	unsigned ptr[4], max = 0, i, j;
	signed ret = 0, score;

	for (i = 0; i < 3; ++i) {
		if (state->bumps[i]) {
			if (!bump_axis_bb (state, balls, i)) return 0;
			if (!state->cnt[i]) goto retn;
			end[max] = state->bumps[i] + state->cnt[i];
			cur[max++] = state->bumps[i];
		}
		ptr[i] = i;
	}
	for (i = 1; i < max; ++i) for (
		j = i; j && cpair_cmp (cur[ptr[j]], cur[ptr[j - 1]]) < 0; --j
	) {
		ptr[3] = ptr[j - 1]; ptr[j - 1] = ptr[j]; ptr[j] = ptr[3];
	}

	for (j = 0;;) {
		struct cpair tmp = *cur[ptr[0]];
		if (!cpair_cmp (&cp, &tmp)) ++j;
		else {
			cp = tmp;
			j = 1;
		}
		// `ret' is only used in narrow_count(), so add 1 instead of `score'.
		if (j == max && (*narrow) (cr, cp.idx[0], cp.idx[1], &score)) ++ret;

		// Also handle duplicates.
		do if (++cur[ptr[0]] == end[ptr[0]]) {
			// See above for why `++ret' is used.
			for (i = 1; i < max && !cpair_cmp (&cp, cur[ptr[i]]); ++i) if (
				++j == max && (*narrow) (cr, cp.idx[0], cp.idx[1], &score)
			) ++ret;
			goto retn;
		} while (!cpair_cmp (&tmp, cur[ptr[0]]));

		for (
			i = 1;
			i < max && cpair_cmp (cur[ptr[i]], cur[ptr[i - 1]]) < 0;
			++i
		) {
			ptr[3] = ptr[i - 1]; ptr[i - 1] = ptr[i]; ptr[i] = ptr[3];
		}
	}

	retn:
	*ret_ = ret;
	return 1;
}
#endif

void cryst_step (struct crystal *cr, float len, unsigned flag) {
	if (!flag) {
		struct cr_var *cell_var = cr->cell->var + cr->vari_perm[cr->nCell[0]];
		cell_var->val[0] += len * cell_var->bound[2];
		if (cell_var->val[0] > cell_var->bound[1]) cell_var->val[0] -= cell_var->bound[2];
		else if (cell_var->val[0] < cell_var->bound[0]) cell_var->val[0] += cell_var->bound[2];
		cell_var->flag = 1;

		if (++cr->nCell[1] == cr->nCell[0]) {
			rng_shuf (cr->rng, cr->vari_perm, cr->nCell[0]);
			cr->nCell[1] = 0;
		}
	}
	else if (flag==1) {
		unsigned n = cr->perm[cr->nPos[1]];
		if (n < cr->nPos[0]) {
			struct cr_pos *pos = cr->poss + n;
			pos->pos[0] += len * pos->bound[2];
			if (pos->pos[0] > pos->bound[1]) pos->pos[0] -= pos->bound[2];
			else if (pos->pos[0] < pos->bound[0]) pos->pos[0] += pos->bound[2];
			pos->flag = 1;
		}
		else if (n > cr->nPos[0]) {
			struct cr_rot *rot = cr->rots + n - cr->nPos[0]; 
			rot->rot[0] += len * rot->bound[2];
			if (rot->rot[0] > rot->bound[1]) rot->rot[0] -= rot->bound[2];
			else if (rot->rot[0] < rot->bound[0]) rot->rot[0] += rot->bound[2];
			rot->flag = 1;
		}
		if (++cr->nPos[1] == cr->nPos[0]) {
			rng_shuf (cr->rng, cr->perm, cr->nPos[0]);
			cr->nPos[1] = 0;
		}
	}
}
/*
void cryst_step_cell (struct crystal *cr, float len) {
	struct cr_var *cell_var = cr->cell->var + cr->vari_perm[cr->nCell[0]];
	cell_var->val[0] += len * cell_var->bound[2];
	if (cell_var->val[0] > cell_var->bound[1]) cell_var->val[0] -= cell_var->bound[2];
	else if (cell_var->val[0] < cell_var->bound[0]) cell_var->val[0] += cell_var->bound[2];
	cell_var->flag = 1;

	if (++cr->nCell[1] == cr->nCell[0]) {
		rng_shuf (cr->rng, cr->vari_perm, cr->nCell[0]);
		cr->nCell[1] = 0;
	}
}
*/

void cryst_ack (struct crystal *cr, int keep) {
	struct cr_rot *rots = cr->rots;
	struct cr_pos *poss = cr->poss;
	struct cr_var *cell_var = cr->cell->var;
	for (unsigned i = 0; i < cr->nPos[0]; ++i) if (poss[i].flag) {
		if (keep) {
			poss[i].pos[1] = poss[i].pos[0];
			poss[i].flag = 0;
		} else poss[i].pos[0] = poss[i].pos[1];
	}
	for (unsigned i = 0; i < 3*cr->nPoly[1]; ++i) if (rots[i].flag) {
		if (keep) {
			rots[i].rot[1] = rots[i].rot[0];
			rots[i].flag = 0;
		} else rots[i].rot[0] = rots[i].rot[1];
	}
	for (unsigned i = 0; i < cr->cell->n; ++i) if (cell_var->flag) {
		if (keep) {
			cell_var[i].val[1] = cell_var[i].val[0];
			cell_var[i].flag = 0;
		} else cell_var[i].val[0] = cell_var[i].val[1];
	}


}

static void theta_eval (struct crystal *cr){
	float t_1,V,a = cr->cell->para[0]; float b = cr->cell->para[1]; float c = cr->cell->para[2];
	float alpha = cr->cell->para[3]*Pi/180; float beta = cr->cell->para[4]*Pi/180; float gama = cr->cell->para[5]*Pi/180;
	V=a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gama)*cos(gama)+2*cos(alpha)*cos(beta)*cos(gama));
	float h,k,l,q;
	int lower_bound = cr->origin/cr->inter;
	struct cryst_surface *surf = cr->surf;
	for (unsigned i = 0; i < cr->nSurf ; ++i) {
		h = surf->hkl[0]; k = surf->hkl[1]; l = surf->hkl[2];
		t_1 = ((h*h)*(b*b)*(c*c)*(sin(alpha)*sin(alpha))
                +(k*k)*(a*a)*(c*c)*(sin(beta)*sin(beta))
                +(l*l)*(a*a)*(b*b)*(sin(gama)*sin(gama))
                +2*h*k*a*b*(c*c)*(cos(alpha)*cos(beta)-cos(gama))
                +2*k*l*(a*a)*b*c*(cos(beta)*cos(gama)-cos(alpha))
                +2*h*l*a*(b*b)*c*(cos(alpha)*cos(gama)-cos(beta)));
		q = sqrt(t_1/V/V);
		if (q*1.54/2<=1) {surf->theta = asin(q*1.54/2)*360/Pi;}
		else surf->theta = cr->origin - cr->inter;
		surf->pos->n = (int)(surf->theta/cr->inter-lower_bound);
		if (surf->pos->n < cr->num && surf->pos->n >= 0){
			surf->flag = 1;
			for (unsigned j = 0; j < cr->nElem; ++j){
				float *asfs = cr->asfs + j*cr->nSurf + i;
				float *para = cr->asfs_para + 9*j;
				*asfs = para[8] + para[0]*exp(-para[4]*q*q/4) 
					+ para[1]*exp(-para[5]*q*q/4)
					+ para[2]*exp(-para[6]*q*q/4)
					+ para[3]*exp(-para[7]*q*q/4);
				
			}
		}
		else {
			surf->flag = 0;
		}
		++surf;
	}
}

static void ful_spec_eval(struct crystal *cr) {
	struct cryst_surface *surf = cr->surf;
	float *cursor = cr->cal_height;
	unsigned i,j;
	//memset function comes from string.h, which can be replaced by loop
	memset(cr->cal_height,0.0,cr->num*sizeof(float));
	for (i = 0; i < cr->nSurf; ++i) if(surf->flag) {
		cursor = cr->cal_height + (surf->pos->n-(Val_spec/2)); // I need add int
		float sc[2] = {
				(float) surf->sc[0] / cr->scales[0],
				(float) surf->sc[1] / cr->scales[0]
			};
		float I = (sc[0] * sc[0] + sc[1] * sc[1]) * surf->mutil;
		for (j = 0; j < Val_spec; ++j){
			*cursor += surf->pos->proj[j] * I ;
			++cursor;
		}
		surf++;
	}
/*
	float max = 0;
	cursor = cr->cal_height;
	for (i = 0; i < cr->num; ++i){
		if (max < *cursor){
			max = *cursor;
		}
		++cursor;
	}

	cursor = cr->cal_height;
	for (i = 0; i < cr->num; ++i){
		*cursor /= max;
		++cursor;
	}
*/
	float max = 0.01;
	for (i = 0; i < cr->num; ++i){
		if (max < cr->cal_height[i]) max = cr->cal_height[i];
	}

	for (i = 0; i < cr->num; ++i){
		cr->cal_height[i] /= max;
	}
}
/*
static void spec_eval (struct crystal *cr) {
	for (unsigned i = 0; i < cr->nHill[0]; ++i) {
		struct cr_hill *hill = cr->hills + i;
		hill->val = 0.0;
		for (unsigned j = 0; j < hill->hill.n; ++j) {
			float sc[2] = {
				(float) hill->peaks[j].sc[0] / cr->scales[0],
				(float) hill->peaks[j].sc[1] / cr->scales[0]
			};
			hill->val += hill->peaks[j].peak.weight *
				(sc[0] * sc[0] + sc[1] * sc[1]);
		}
	}
}
*/
#ifdef BENCHMARK
static void spec_eval_bc (struct crystal *cr) {
	for (unsigned i = 0; i < cr->nHill[0]; ++i) {
		struct cr_hill *hill = cr->hills + i;
		hill->val = 0.0;
		for (unsigned j = 0; j < hill->hill.n; ++j) {
			float *sc = hill->peaks[j].bsc;
			hill->val += hill->peaks[j].peak.weight *
				(sc[0] * sc[0] + sc[1] * sc[1]);
		}
	}
}
#endif

// use the definition of R_bragg, R_w(not defined)
//R_b, the remain factor < 0.1
/*
static float Rb_eval (struct crystal *cr) {
	float ret_num = 0.0, ret_den = 0.0, I_obs, I_cal; // ret_num is the Numerator, ret_den is denominator
	unsigned i;
	for (i = 0; i< cr->num; ++i){
		I_cal = cr->cal_height[i];
		I_obs = cr->height[i];
		ret_num += (I_obs-I_cal)*(I_obs-I_cal);
		ret_den += I_obs*I_obs;
	}
	return sqrt(ret_num/ret_den);
}
*/
static float Rw_eval (struct crystal *cr) {
	float ret_num = 0.0, ret_den = 0.0, I_obs, I_cal, limit=0; // ret_num is the Numerator, ret_den is denominator
        unsigned i;
        for (i = 0; i< cr->num; ++i) {
		if (limit > cr->height[i]) limit=cr->height[i];
	}
        for (i = 0; i< cr->num; ++i){
                I_cal = cr->cal_height[i];
                I_obs = cr->height[i];
                //if (I_obs+I_cal<0.06) continue;
                ret_num += (I_obs-I_cal)*(I_obs-I_cal);//(I_obs-limit+0.01);
                ret_den += I_obs*I_obs;//(I_obs-limit+0.01);
        }
        return sqrt(ret_num/ret_den);
}

static float metric_dist (
	metric const *met, float const a[3], float const b[3]
) {
	float dif[3];
	for (unsigned i = 0; i < 3; ++i) dif[i] = b[i] - a[i];
	return sqrt (metric_get (met, dif));
}

#ifdef BENCHMARK
int narrow_orig (
	struct crystal const *cr, unsigned i, unsigned j, signed *ret_
) {
	struct cr_ball const *balls = cr->atoms->scats->balls;
	float ret = metric_dist (cr->met, balls[i].xyz, balls[j].xyz) /
		balls[i].dists[balls[j].elem];
	*ret_ = cr->scales[1] * (ret > cr->ctl[1] ? 0.0 : (
			ret < cr->ctl[0] ? 1.0 : cr->ctl[2] * (cr->ctl[1] - ret)
		));
	return ret < cr->ctl[1];
}

#else
static
#endif
int narrow_cryst (
	struct crystal const *cr, unsigned i, unsigned j, signed *ret_
) {
	struct cr_ball const *balls = cr->atoms->scats->balls;
	float ret = metric_dist (cr->met, balls[i].xyz, balls[j].xyz) /
		balls[i].dists[balls[j].elem];
	*ret_ = cr->scales[1] * (balls[i].n + balls[j].n) *
		(ret > cr->ctl[1] ? 0.0 : (ret < cr->ctl[0] ? 1.0 :
			cr->ctl[2] * (cr->ctl[1] - ret)
		));
	return ret < cr->ctl[1];
}

static int bump_eval (struct crystal *cr, narrow_f *narrow) {
	struct rbforest *bump = &cr->bump;
	struct cr_ball *balls = cr->atoms->scats->balls;
	signed score;
	unsigned i, j;
	for (i = 0; i < bump->n; ++i) {
		if (balls[i].flag & 0x1) rbforest_rm (bump, i);
	}
	for (i = 0; i < bump->n; ++i) {
		unsigned flag0 = balls[i].flag;
                unsigned ind0 = balls[i].ind;
		if (!(flag0 & 0x2)) continue;
		for (j = 0; j < bump->n; ++j) {
			unsigned flag1 = balls[j].flag;
                        unsigned ind1 = balls[j].ind;
			if (i == j || (flag1 & 0x2 && i > j) || !(
				(flag0 | flag1) & 0x1 && (*narrow) (cr, i, j, &score)
			)) continue;
                        if ((ind1 == ind0) && (ind0 != 0)) {
				continue;
			}
			if (!rbforest_add (bump, i, j, score)) return 0;
		}
	}
	return 1;
}

static void score_eval (struct crystal *cr) {
	cr->scores[1] = Rw_eval(cr);
	if ((cr->scores[2] /= 2 * (cr->nAtom[1] + cr->nPoly[3])) > 1.0) cr->scores[2] = 1.0;
	cr->scores[0] = cr->scales[2] * cr->scores[2] +
		(1 - cr->scales[2]) * cr->scores[1];
}



//lyt add 2021/11/16
static float Gauss(float theta, float H) {
	// 0.9394372787 is sqrt(4Ln2/Pi)
	return 0.9394372787*exp(-4*Ln2*(theta)*(theta)/H/H)/H;
}

static float Los(float theta, float H) {
	return 2/(1+(4*(theta)*(theta)/H/H))/H/Pi;
}
// this is the second cw profile function in GASA, which is pV function
static void cwfun_2(float theta, float *proj) {
	unsigned i;
	float scale = 4;
        float  H_g, H_l, n, L_g, L,GU = 6*scale, GV = -6*scale, GW = 5*scale, Lx = 4*scale, Ly = 4*scale;
        float U = GU/1803.4, V = GV/1803.4, W = GW/1803.4, X = Lx/100, Y = Ly/100;
        float tan_t = tan(theta*Pi/360), cos_t = cos(theta*Pi/360);
        H_g = sqrt(U*tan_t*tan_t+V*tan_t+W);
        H_l = (X/cos_t) + (Y)*tan_t;
        L_g = sqrt(8*Ln2)*H_g;
	/*
	float U = 4, V = -4, W = 5, Lx = 4, Ly = 4, H_g, H_l, n, L_g, L;
	float tan_t = tan(theta*Pi/360), cos_t = cos(theta*Pi/360);
	H_g = sqrt((U*tan_t*tan_t+V*tan_t+W)/cos_t/cos_t);
	H_l = (Lx/cos_t) + (Ly)*tan_t;
	L_g = sqrt(8*Ln2)*H_g;  // this is FWHM = ..*sigma； FWHM = 2*sqrt(ln(2))*sigma
	*/
	L = pow(pow(L_g,5)+2.69269*pow(L_g,4)*H_l \
			+2.42843*pow(L_g,3)*H_l*H_l +4.47163*pow(L_g,2)*pow(H_l,3) \
			+0.07842*pow(L_g,1)*pow(H_l,4)+pow(H_l,5),0.2);
	
	n = 1.36603*(H_l/L) - 0.47719*pow(H_l/L,2) + 0.11116*pow(H_l/L,3);
	
	float theta_re_0 = 0.02*(50*theta-(int)(50*theta)); //the intern is 0.02
	
	float theta_re = -theta_re_0 - 0.02*(int)(Val_spec/2) - 0.02; // I add minus theta_re_0
	for (i=0; i<Val_spec; ++i){
		theta_re += 0.02;
		proj[i] = n*Los(theta_re,H_l)+(1-n)*Gauss(theta_re,H_g);
	}
/*
	float a[10] = {Gauss(-0.1,H_g),Gauss(-0.08,H_g),Gauss(-0.06,H_g), \
					Gauss(-0.04,H_g),Gauss(-0.02,H_g),Gauss(-0.0,H_g), \
					Gauss(0.02,H_g),Gauss(0.04,H_g),Gauss(0.06,H_g), \
					Gauss(0.08,H_g),Gauss(0.10,H_g)};
	float b[10] = {Los(-0.1,H_l),Los(-0.08,H_l),Los(-0.06,H_l), \
					Los(-0.04,H_l),Los(-0.02,H_l),Los(-0.0,H_l), \
					Los(0.02,H_l),Los(0.04,H_l),Los(0.06,H_l), \
					Los(0.08,H_l),Los(0.10,H_l)};
*/
}

//lyt add 2021/11/6
int peak_fun (crystal *cr) {
	unsigned i;
	struct cryst_surface *hkl_val = cr->surf;
	for (i = 0; i < cr->nSurf; ++i) if (hkl_val->flag){
		if (i) hkl_val++; //hkl_val->pos->n = (int)(hkl_val->theta*50)-1;
		cwfun_2(hkl_val->theta,hkl_val->pos->proj);
	}

	if (!(hkl_val->pos->proj)) goto profile_err;
	return 1;
	free(hkl_val->pos->proj);profile_err: 
	return 0;
}

int cryst_eval_0 (struct crystal *cr) {
	struct cr_poly *poly = cr->polys;
	struct cr_rot *rot = cr->rots;
	struct cr_pos *pos = cr->poss;
	struct cr_atom *atom = cr->atoms;
	struct cr_ball *ball = atom->scats->balls;
	struct cr_var *cell_var = cr->cell->var;
	struct linear *cell_lin = cr->cell->lin;
	float *angles = cr->cell->para+3, *abc = cr->cell->para; 
	float new_para;
	unsigned i;
	for (i = 0; i < cr->cell->n; ++i) if (cell_var->flag){
		for (unsigned j = 0; j < 6; ++j){
			new_para = cell_var->val[0]*cell_lin->axis + cr->cell->lin[j].zoom;
			cr->cell->para[j] = new_para>0 ? new_para : cr->cell->para[j];
			cell_var->flag = 0;
			++cell_lin; 
		}
		++cell_var; 
		cr->cell->flag = 1;
	}
	if (cr->cell->flag){
		theta_eval (cr);
		peak_fun (cr);
		cr->cell->flag = 0;
	}
	for (i = 0; i < cr->nPoly[1]*3; ++i, ++rot) if (rot->flag) {
		rot->poly->angles[rot->axis] = torus_pos ((rot->rot[0]-rot->bound[0])/rot->bound[2]) * rot->bound[2] + rot->bound[0];
		rot->poly->flag = 1;	
	}

	for (i = 0; i < cr->nPos[0]; ++i, ++pos) if (pos->flag) {
		pos->atom->xyz[pos->axis] = torus_pos (pos->pos[0]);
		pos->atom->flag = 1;
	}
        for (i = 0; i < cr->nPoly[1]; ++i, ++poly) if (poly->flag || poly->orig->flag) {
		struct cr_atom *extra_atom = poly->atom;
		if (poly->flag) rot_mat_mk(poly->trans, abc, angles, poly->angles);
		for (unsigned j = 0; j < poly->n; ++j, ++extra_atom) {
			mat_mult_add(poly->trans, poly->ref_coor[j], poly->orig->xyz, extra_atom->xyz);
			extra_atom->flag = 1;
		}
		poly->flag = 0;	
	}
	for (i = 0; i < cr->nAtom[0]+cr->nPoly[2]; ++i, ++atom) {
		surf_eval (cr, atom);
		/*
                int cos_kx=0,sin_kx=0;
                int k;
                for (k = 0; k < 24; ++k) {
                        cos_kx += (cr->atoms->scats+k)->scs[0][0];
                        sin_kx += (cr->atoms->scats+k)->scs[0][1];
                        printf("%d,%d,this is %d th \n",cos_kx,sin_kx,k+1);
                }
                */
		atom->flag = 0;
	}
	
	ful_spec_eval (cr);
	//spec_eval (cr);
	/*
	for (i = 0; i<cr->num; i++)
		printf("%f\n", cr->cal_height[i]);
	*/	
	if (cr->scales[2] != 0.0) {
		if (!bump_eval (cr, &narrow_cryst)) return 0;
		for (i = 0; i < cr->nAtom[1] + cr->nPoly[3]; ++i, ++ball) ball->flag &= 0x2;
	}
	cr->scores[2] = cr->bump.score / cr->scales[1];
	score_eval (cr);
	//printf("this is position and cell_parameter: %f %f %f %f %f %f %f %f %f \n", cr->cell->para[0], cr->cell->para[1], cr->cell->para[2],
                       //cr->cell->para[3], cr->cell->para[4], cr->cell->para[5], cr->atoms->xyz[0], cr->atoms->xyz[1], cr->atoms->xyz[2]);
	return 1; 
}

int cryst_eval (struct crystal *cr) {
	struct cr_poly *poly = cr->polys;
	struct cr_rot *rot = cr->rots;
	struct cr_pos *pos = cr->poss;
	struct cr_atom *atom = cr->atoms;
	struct cr_ball *ball = atom->scats->balls;
	struct cr_var *cell_var = cr->cell->var;
	struct linear *cell_lin = cr->cell->lin;
	float *angles = cr->cell->para+3, *abc = cr->cell->para; 
	float new_para;
	unsigned i;
	for (i = 0; i < cr->cell->n; ++i) if (cell_var->flag){
		for (unsigned j = 0; j < 6; ++j){
			new_para = cell_var->val[0]*cell_lin->axis + cr->cell->lin[j].zoom;
			cr->cell->para[j] = new_para>0 ? new_para : cr->cell->para[j];
			cell_var->flag = 0;
			++cell_lin; 
		}
		++cell_var;
		cr->cell->flag = 1;
	}
	if (cr->cell->flag){
		theta_eval (cr);
		peak_fun (cr);
		cr->cell->flag = 0;
	}

	//0.26 0.107932 0.408015 0.5422 1.01065 0.6615
/*	(pos)->flag = 1;	pos->pos[0] = 0.26;
	(pos+1)->flag = 1;	(pos+1)->pos[0] = 0.107932;
	(pos+2)->flag = 1;	(pos+2)->pos[0] = 0.40815;
	(pos+3)->flag = 1;	(pos+3)->pos[0] = 0.6078;
	(pos+4)->flag = 1;	(pos+4)->pos[0] = 0.3922;
	(pos+5)->flag = 1;	(pos+5)->pos[0] = 0.2898;
	(pos+6)->flag = 1;	(pos+6)->pos[0] = 0.439622;
	(pos+7)->flag = 1;	(pos+7)->pos[0] = 0.204138;
	(pos+8)->flag = 1;	(pos+8)->pos[0] = 0.103251;
	//(pos+9)->flag = 1;	(pos+9)->pos[0] = 0.132471;
	rot->flag = 1; rot->rot[0] = 80.3155;
	(rot+1)->flag = 1; (rot+1)->rot[0] = 0.00579321;
	(rot+2)->flag = 1; (rot+2)->rot[0] = 119.208;
	(rot+3)->flag = 1; (rot+3)->rot[0] = 59.9813;
	(rot+4)->flag = 1; (rot+4)->rot[0] = 179.961;
	(rot+5)->flag = 1; (rot+5)->rot[0] = 28.3855;
*/	
	for (i = 0; i < cr->nPos[0]; ++i, ++pos) if (pos->flag) {
		pos->atom->xyz[pos->axis] = torus_pos (pos->pos[0]);
		pos->atom->flag = 1;
	}
	for (i = 0; i < cr->nPoly[1]*3; ++i, ++rot) if (rot->flag) {
		rot->poly->angles[rot->axis] = torus_pos ((rot->rot[0]-rot->bound[0])/rot->bound[2]) * rot->bound[2] + rot->bound[0];
		rot->poly->flag = 1;	
	}

	for (i = 0; i < cr->nPoly[1]; ++i, ++poly) if (poly->flag || poly->orig->flag) {
		if (poly->orig->flag) {	
			unsigned nn[2];
			balls_eval (cr, poly->orig, nn);
			}
                struct cr_atom *extra_atom = poly->atom;
		if (poly->flag) rot_mat_mk(poly->trans, abc, angles, poly->angles);
		for (unsigned j = 0; j < poly->n; ++j, ++extra_atom) {
			mat_mult_add(poly->trans, poly->ref_coor[j], poly->orig->scats->balls->xyz, extra_atom->xyz);
			extra_atom->flag = 1;
		}
		poly->flag = 0;	
	}
	for (i = 0; i < cr->nAtom[0]+cr->nPoly[2]; ++i, ++atom) if (atom->flag) {
		surf_eval (cr, atom);
		atom->flag = 0;
	}
	ful_spec_eval (cr);
	
	//for (i = 0; i<cr->num; i++)
	//	printf("%f\n", cr->cal_height[i]);
	
	if (cr->scales[2] != 0.0) {
		if (!bump_eval (cr, &narrow_cryst)) return 0;
		for (i = 0; i < cr->nAtom[1] + cr->nPoly[3]; ++i, ++ball) 
			{ball->flag &= 0x2;}
	}
	cr->scores[2] = cr->bump.score / cr->scales[1];
	score_eval (cr);
	return 1; 
}

#ifdef BENCHMARK
static void pre_eval (struct crystal *cr, int sc) {
	struct cr_peak *peak = cr->hills->peaks;
	struct cr_pos *pos = cr->poss;
	unsigned i;
	for (i = 0; i < cr->nPos[0]; ++i, ++pos) {
		pos->atom->xyz[pos->axis] = pos->pos[0];
	}
	for (i = 0; i < cr->nHill[1]; ++i, ++peak) memset
		(peak->bsc, 0, 2 * sizeof (float));
	for (i = 0; i < cr->nAtom[0]; ++i) atom_eval_bb (cr, cr->atoms + i, sc);
}

int cryst_eval_bb (struct crystal *cr, broad_f *broad, narrow_f *narrow) {
	pre_eval (cr, 0);
	if (broad && narrow) {
		signed score;
		if (!(*broad) (cr, narrow, &score)) return 0;
		cr->scores[0] = score;
	}
	return 1;
}

int cryst_eval_bc (struct crystal *cr) {
	signed score;
	pre_eval (cr, 1);
	spec_eval_bc (cr);
	if (!bump_eval_bc (cr, &narrow_cryst, &score)) return 0;
	cr->scores[2] = score / cr->scales[1];
	score_eval (cr);
	return 1;
}
#endif

unsigned cryst_dof (struct crystal const *cr) {
	return cr->nPos[0] + cr->nPoly[1]*3; //+ cr->cell->n;
}

void cryst_dump (struct crystal const *cr, uint32_t buf[]) {
	unsigned i;
	/*
	for (i = 0; i < cr->nCell[0]; ++i) {
		*buf = hftonl (cr->cell->var[i].val[0]);
		buf++;
	}
	*/
	for (i = 0; i < cr->nPos[0]; ++i) {
		*buf = hftonl (cr->poss[i].pos[0]);
		buf++;
	}
        for (i = 0; i < cr->nPoly[1]*3; ++i) {
		*buf = hftonl (cr->rots[i].rot[0]);
		buf++;
	}
}

int cryst_load (struct crystal *cr, uint32_t const buf[]) {
	struct cr_pos *poss = cr->poss;
	for (unsigned i = 0; i < cr->nPos[0]; ++i) {
		poss[i].pos[0] = nltohf (buf[i]);
		poss[i].flag = 1;
		if (
			!isfinite (poss[i].pos[0]) ||
			poss[i].pos[0] < poss[i].bound[0] ||
			poss[i].pos[0] > poss[i].bound[1]
		) return 0;
	}
	return 1;
}

float const *cryst_scores (struct crystal const *cr) {
	return cr->scores;
}

int mc_grid (struct crystal *cr) {
        float i,j;
        for (i=0;i<1;i+=0.001) {
                for (j=0;j<1;j+=0.001){
                                cr->atoms->xyz[0] = i;
                                cr->atoms->xyz[1] = i;
                                cr->atoms->xyz[2] = j;
                                cryst_eval_0(cr);
                                printf("%f %f %f \n",i,j,cr->scores[1]);
                }
        }
        return 1;
}
//extern in optim.h
int mc_bias (struct crystal *cr) {
	cryst_eval_0(cr);
	fprintf(stdout,"%f",cr->scores[1]);
	return 1;
}

//extern in cryst.h
int init_scan_var(crystal *cr, int argc, char const *const *argv){
	struct cr_pos *pos = cr->poss;
	struct cr_var *cell_var = cr->cell->var;
	if (argc!=(cr->nPos[0]+cr->cell->n+1)) goto num_err;
	for (int i=0; i<cr->cell->n; i++){
		if (!scan2_float (argv[i+1], &cell_var->val[0])) {
		fprintf (stderr, "Usage: decr_mcs n cell_parameters n coordinates < cryst.cr\n");
		goto num_err;}
		cell_var->flag = 1;
		++cell_var;
	}
	for (int i=0;i<cr->nPos[0];i++){
		if (!scan2_float (argv[i+1+cr->cell->n], &pos->pos[0])) {
		fprintf (stderr, "Usage: decr_mcs n cell_parameters n coordinates < cryst.cr\n");
		goto num_err;
		}
		pos->flag = 1;
		pos++;
	}
	return 1;num_err: 
	return 0;
}

#ifdef BENCHMARK
unsigned cryst_nball (struct crystal const *cr) {
	return cr->nAtom[1];
}

void cryst_set (struct crystal *cr, float const xyzs[][3]) {
	struct cr_pos *poss = cr->poss;
	for (unsigned i = 0; i < cr->nPos[0]; ++i) {
		poss[i].pos[0] = xyzs[poss[i].atom - cr->atoms][poss[i].axis];
		poss[i].flag = 1;
	}
}

void cryst_write (struct crystal const *cr) {
	printf ("%g %g\n", cr->scores[1], cr->scores[2]);
	for (unsigned i = 0;;) {
		float *xyz = cr->atoms[i].scats->balls->xyz;
		printf ("%g,%g,%g", xyz[0], xyz[1], xyz[2]);
		if (++i < cr->nAtom[0]) printf (" ");
		else {
			printf ("\n");
			break;
		}
	}
}
#endif
