typedef struct metric metric;

// Inputs are (a, b, c) and (alpha, beta, gamma).  The angles shall be in
// degrees; the cell shall be chosen so that the angles are as close to 90 as
// possible (as crystallographers usually do).  Returns NULL iff failed;
// otherwise returned pointer shall be free()d after use.
extern metric *metric_mk (float const abc[3], float const angles[3]);

// Inputs are (a, b, c) and (alpha, beta, gamma) and angles of rotations
extern void rot_mat_mk (float trans[3][3], float const abc[3], float const angles[3], float dcm_angles[3]);

// compute the transformation of coordinates
extern void mat_mult_add (float const m[3][3], float const x[3], float const y[3], float ans[3]);

// Computes d^2 from torus displacement (x, y, z) in (-1, 1)^3.
extern float metric_get_slow (metric const *met, float const xyz[3]);
// The fast variant is only guaranteed to be correct if
// sum (abs (x - round (x)) for x in xyz) < 0.5.
extern float metric_get_fast (metric const *met, float const xyz[3]);
// Automatically switches to the slow variant only if necessary, useful when
// the slow variant is rarely (but still possibly) needed, which is often the
// case.
extern float metric_get (metric const *met, float const xyz[3]);

