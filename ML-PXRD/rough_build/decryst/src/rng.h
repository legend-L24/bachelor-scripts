#include <stdint.h>

typedef struct rng_t rng_t;
typedef struct step_t step_t;

// Returns NULL iff failed; otherwise returned pointer shall be free()d after
// use.
extern rng_t *rng_mk (uint64_t seed);
// Uniform in 0 : 2^32.
extern uint32_t rng_next (rng_t *rng);
// Uniform in [0, 1).
extern double rng_float (rng_t *rng);
// Uniform in 0 : n (> 0).
extern unsigned rng_dice (rng_t *rng, unsigned n);
// Biased in 0 : n (> 0) according to upper quantiles (ending with 1) in `qs'.
extern unsigned rng_bdice (rng_t *rng, float const qs[], unsigned n);
// Shuffles an `arr'ay of length `n' (> 0).
extern void rng_shuf (rng_t *rng, unsigned arr[], unsigned n);

// The step length distribution; initialise with step_mag() before use.
extern step_t *step_mk (rng_t *rng);
// Initialises or modifies the `mag'nitude constant.
extern void step_mag (step_t *step, float mag);
// A random variable from the distribution with range (-0.5, 0.5).
extern float step_get (step_t const *step);

// Initialises `seed', should be MT-safe; returns 1 on success, 0 on failure.
extern int seed_mk (uint64_t *seed);
// Combination of seed_mk() and rng_mk().
extern rng_t *rng_mk2 (void);
// Converts weights (>= 0, and not all 0) in `qs' to upper quantiles.
extern void wstoqs (float qs[], unsigned n);

