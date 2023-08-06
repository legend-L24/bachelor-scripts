#include <stdint.h>

// J. Lam, J.-M. Delosme. "An Efficient Simulated Annealing Schedule:
// Implementation and Evaluation". Technical Report 8817, Department of
// Electrical Engineering, Yale University, 1988.
// K.-W. Chu, Y.-F. Deng, J. Reinitz. "Parallel Simulated Annealing by Mixing of
// States". Journal of Computational Physics, 1999, 148(2): 646-662.

typedef struct sa_sched sa_sched;
typedef struct sa_comp sa_comp;

// Gathers statistics of `n' (> 1) states (including the current state, assuming
// a auccessful cryst_eval() / cryst_ack() pair beforehand), setting `ret' to
// (mean, sd, best, tv_of_best, cd_of_best); returns 1 on success, 0 on failure.
extern int mc_stat (crystal *cr, rng_t *rng, unsigned n, float ret[5]);

// use grid search to generate R value
extern int mc_grid (crystal *cr);

// try to evaluate one coordinate
extern int mc_bias(crystal *cr);

// Initialises control parameters based on `ab' (alpha, beta) and `init'
// (lambda, mean, sd), assuming s_0 = 0; sets `ctl' to (a, b, d, e, step
// magnitude, r, s) with (mag, r, s) hardcoded to (1, 0, 0); returns NULL iff
// failed; otherwise returned pointer shall be free()d after use.
// r: 4 * lambda * rho_0 * ((1 - rho_0) / (2 - rho_0)) ^ 2.
extern sa_sched *sa_sched_mk
	(float const ab[2], float const init[3], float ctl[7]);
// Updates the parameters based on `stat' (mean, sd, latest s, rho_0); returns 1
// on success, or 0 if unhandled Inf/NaN occured.
extern int sa_sched_step (sa_sched *sched, float const stat[4], float ret[6]);

// Returns NULL iff failed.
extern sa_comp *sa_comp_mk (crystal *cr, rng_t *rng);
extern void sa_comp_fin (sa_comp *comp);
// Scores (nc, tv, cd) and degrees of freedom for the best model known.  The
// arrays will be valid as long as `comp' is valid, and reflect value changes.
extern float const *sa_comp_best (sa_comp *comp);
extern uint32_t const *sa_comp_bbuf (sa_comp *comp);

// Runs tau = n[0] (> 0) moves; assuming a successful cryst_eval() / cryst_ack()
// pair beforehand, and that this is the n[1]-th (0-based) of n[2] (> 0)
// parallel processes; with initial parameters `ctl' (a, b, d, e, mag, r, latest
// s after moves on all cores, delta) (s >= 0); sets `iret' to (num of accepted
// moves), `fret' to (score sum, sum of squared deviations wrt estimated means,
// final inverse temperature); returns 1 on success, 0 on failure.
// delta (> 1): s[i+1] (counting only moves on this core) will be set to
// delta * s[i] if the latter is exceeded.
// if flag = 1, otpimizer the coordinates of atoms; if flag = 0, optimizer the parameters of cell
extern int sa_comp_get (
	sa_comp *comp, unsigned const n[3], float const ctl[8],
	unsigned iret[1], float fret[3], unsigned flag
);

