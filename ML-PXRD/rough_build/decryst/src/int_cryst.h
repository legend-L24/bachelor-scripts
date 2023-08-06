struct cr_peak {
	struct peak_t peak;
#ifdef BENCHMARK
	float bsc[2];
#endif
	signed sc[2];
};

struct cr_hill {
	struct hill_t hill;
	struct cr_peak *peaks;
	float val;
};
// I need add some new part

struct cr_ball {
	float *dists, xyz[3];
#ifdef BENCHMARK
	float bound[3][2];
#endif
	unsigned flag, elem, n, ind;
};

struct cr_scat {
	struct cr_ball *balls;
	signed (*scs)[2];
};
struct cr_atom {
	struct affine *affs;
	struct cr_scat *scats;
	float *asfs, xyz[3];
#ifdef BENCHMARK
	float *wid;
#endif
	unsigned (*aids)[3], wyck[3];
	int flag;
};

struct cr_poly {
	struct cr_atom *orig;
	struct cr_atom *atom;
	float angles[3], (*ref_coor)[3];
        float trans[3][3];
	int flag;
	unsigned n;
	// which is DCM rotation matrix 
};


struct cr_rot {
	float rot[2], bound[3];
	struct cr_poly *poly;
	unsigned axis, flag;
};

struct cr_pos {
	float pos[2], bound[3];
	struct cr_atom *atom;
	unsigned axis, flag;
};

struct cr_var {
	float val[2], bound[3];
	unsigned flag;
};


struct cr_cell {
	float *para; // this is a,b,c, alpha,beta,gama
	struct cr_var *var;
	unsigned n;
	int flag;
	struct linear *lin; // the lin[n*6] n n n ...
};
struct cr_peak_para {
	float *para;	//U, V, W, X, Y are five parameters of cw2_fun
	struct cr_var *var;
	unsigned n;
};

//the experimental spectrm is expressed by num, height[num], origin(about theta), 
//the caculated spectrum is expressed by num, cal_height[num],the mutil of surface is multiplicity times Lorentz-polarisation factor(peaks->weights)
//asfs_para data is similar to asfs and the nine parameter of atoms scatter factor
struct  crystal {
	struct rbforest bump;
#ifdef BENCHMARK
	struct bump_t *bbump;
#endif
	rng_t *rng;
	metric *met;
	struct cr_peak_para *peak_para;
	struct cr_cell *cell;
        struct cr_pos *poss;
	struct cr_rot *rots;
	struct cr_atom *atoms;
	struct cr_poly *polys;
	struct cr_hill *hills;
	struct affine *affs;
	struct cryst_surface *surf;
	unsigned *perm, *vari_perm, (*aids)[3], nHill[2], nAtom[2], nPoly[4], nPos[2], nCell[2], nOrig, num, nSurf, nElem;
#ifdef BENCHMARK
	float (*wids)[3];
#endif
	float *dists, *asfs, *asfs_para, (*origs)[3], ctl[3], scales[3], scores[3], *height, *cal_height, origin, inter;  
};

#ifndef BENCHMARK
typedef int narrow_f (crystal const *, unsigned, unsigned, signed *);
#endif

extern float torus_pos (float x);
// this is evalute the intensity of every crystal surface
extern void balls_eval (struct crystal *cr, struct cr_atom *atom, unsigned nn[2]);
extern void surf_eval (struct crystal *cr, struct cr_atom *atom);
extern void atom_eval (crystal *cr, struct cr_atom *atom);
#ifdef BENCHMARK
extern void atom_eval_bb (crystal *cr, struct cr_atom *atom, int sc);
#endif
