#include <stdint.h>
#include <stdio.h>
#define Val_spec 35 //range of num is from 10-20
// this struct is added by lyt in 2021 10 31
struct cryst_surface_t {
	signed hkl[3];
	float mutil;
	float theta;
};
struct peak_shape {
	int n;
	float proj[Val_spec];
};

struct cryst_surface {
	signed hkl[3];
	float mutil;
	float theta;
	signed sc[2];
        int flag;
	struct peak_shape *pos;
};

// this struct is added by lyt
struct spec_t {
	float theta;
	float I;
};

struct peak_t {
	signed hkl[3];
	float weight;  // Shall be > 0.
};

struct hill_t {
	unsigned n;  // Number of affiliated peaks, shall be > 0.
	float val0;  // Shall be >= 0.
};

struct linear {
	unsigned axis;
	float zoom;
};

struct affine {
	struct linear lin[2];
	float move;
	unsigned n;
};
struct extra_atom_t {
	unsigned elem;
	float xyz[3];
};
struct atom_t {
	unsigned elem, wyck;
	// Each lower bound shall not be larger than the corresponding upper bound.
	// Set the the bounds on one axis exactly (on the binary representation
	// level) equal to declare the coordinate on that axis as unused or (if
	// used) constant.
	float bound[3][2];
};
struct poly_t {
	unsigned center, wyck, n;
	// Each lower bound shall not be larger than the corresponding upper bound.
	// set the rotation matrix as DCM matrix (in most cases, three angles needs to be
	// optimized
	float bound[3][2];
};
// lyt add the spec_t,nSpec in 2021/10/9
struct crin_t {
	struct peak_t *peaks;
	struct hill_t *hills;
	struct atom_t *atoms;
	struct affine *affs;
        struct spec_t *spec;
        struct poly_t *polys;
        struct cryst_surface_t *surf;
        struct extra_atom_t *extra_atoms;
	float (*origs)[3], *dists, *asfs, abc[3], angles[3], ctl[3], *asfs_para;
	// `wycks': array of (number of affs, number of equivalent positions or
	// centrosymmetric pairs, centrosymmetry of the equivalent positions wrt the
	// origin); nAff and nPos shall be > 0; centrosymmetry shall be <= 1.
	// If a Wyckoff position is centrosymmetric wrt the origin, then only one
	// position would be needed from each centrosymmetric pair.
	unsigned (*wycks)[3], (*aids)[3], nHill, nOrig, nWyck, nElem, nAtom, nSpec, nSurf, sg_number, nPoly;
};

typedef struct crystal crystal;

// Returns NULL iff failed; otherwise all variable coordinates are randomised.
// The collision detection score function f(s) wrt s = d(i, j) / dists[i][j] is
// a piecewise linear function with control nodes (0 <= ctl[1] < ctl[2])
//   f(0) = f(ctl[1]) = 1, f(ctl[2]) = f(+inf) = 0.
// The net cost function is (0 <= ctl[0] < 1)
//   nc = (1 - ctl[0]) * tv + ctl[0] * cd,
// where tv (in [0, 1]) is the total variation distance between computed and
// actual spectra (both based on hills, not peaks), and cd (with a hardcoded
// upper limit at 1) is the sum of collision detection scores (each pair counted
// once) divided by the number of atoms in the whole cell (i.e. not only those
// in the asymmetric unit).
extern crystal *cryst_mk (
	// `peaks' shall be aligned according to `hills', eg. if hills[].n is
	// [2, 4, ...], then peaks[0 : 2] will be affiliated to hills[0],
	// peaks[2 : 6] to hills[1], etc; origs[0] shall be (0, 0, 0).
	struct peak_t const peaks[], struct hill_t const hills[], struct spec_t const spec[],
	float const origs[][3], struct affine const affs[],
	// `affs'/`aids' vs. `wycks' is similar to `peaks' vs. `hills', just note
	// that each Wyckoff letter corresponds to an array of `affine's and an
	// array of `affine'-id 3-tuples.  Each Wyckoff letter's `affine's shall be
	// numbered separately, always 0-based; first equivalent position of each
	// atom shall be (x, y, z) of the atom itself.
	unsigned const aids[][3], unsigned const wycks[][3],
	// `dists' and `asfs' shall be of sizes [nElem][nElem] and [nElem][nPeak],
        //asfs_para is the similar data to asfs,[nElem][nPeak]
	// respectively, and aligned accordingly, with every element > 0;
	// dists[i][j] == dists[j][i] shall exactly hold.
	float const dists[], float const asfs[], float const asfs_para[],struct atom_t const atoms[],
        struct poly_t const polys[], struct extra_atom_t const extra_atoms[], struct cryst_surface_t const surf[],
	float const abc[3], float const angles[3], float const ctl[3],
	// All following `n's shall be > 0.
	rng_t *rng, unsigned nHill, unsigned nOrig,
	unsigned nWyck, unsigned nElem, unsigned nAtom, unsigned nSpec, unsigned nSurf, unsigned sg_number, unsigned nPoly
);
extern void cryst_fin (crystal *cr);

// Reads crystal parameters from open `file'; returns NULL iff failed.
extern struct crin_t *crin_read (FILE *file);
extern void crin_fin (struct crin_t *crin);
// Returns NULL iff failed.
extern crystal *crin_eval (struct crin_t const *crin, rng_t *rng);
// Returns NULL iff failed.
extern crystal *cryst_read (rng_t *rng, FILE *file);

// Allowed modification sequences:
//   mk -> eval, eval -> ack, ack -> step|load, step|load -> eval, !fin -> fin.
// The first call to cryst_ack() after cryst_mk() shall set `keep' to 1.
// cryst_step() shall not be called when the degree of freedom is 0.

// Moves specified step `len'gth (in [-0.5, 0.5]) in position parameters or cell parameter.
// if flag = 1, otpimizer the coordinates of atoms; if flag = 0, optimizer the parameters of cell
extern void cryst_step (crystal *cr, float len, unsigned flag);
// After evaluation, set whether to keep previous move.
extern void cryst_ack (crystal *cr, int keep);
// compute the peak function(only 10-20 valid points); returns 1 on success, 0 on failure.
extern int peak_fun (crystal *cr);
// Computes the score; returns 1 on success, 0 on failure.
extern int cryst_eval (crystal *cr);
// Returns the degree of freedom.
extern unsigned cryst_dof (crystal const *cr);

// Dumps/loads the degrees of freedom from/to a `buf'fer for `float's converted
// to `uint32_t's in the network byte order.  cryst_load() expects each degree
// of freedom to be within the corresponding bound; returns 1 if the condition
// is met, or 0 otherwise.
extern void cryst_dump (crystal const *cr, uint32_t buf[]);
extern int cryst_load (crystal *cr, uint32_t const buf[]);
// The scores (see above for the definitions) of a model: (nc, tv, cd).  The
// array will be valid as long as `cr' is valid, and will reflect value changes.
extern float const *cryst_scores (crystal const *cr);

// in bayes optimization, input x,y,z, output R_eval 
extern int init_scan_var(crystal *cr, int argc, char const *const *argv);
