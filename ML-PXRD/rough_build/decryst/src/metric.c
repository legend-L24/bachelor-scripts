#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "metric.h"
#include "utils.h"
#ifdef BENCHMARK
#include "bench_metric.h"
#endif

struct metric {
	float trans[3][3], vert8[8][3];
#ifdef BENCHMARK
	float g[3][3], sq[8], vert27[27][3];
#endif
};
static void vec_round (float const xyz[3], float r[3]) {
	for (unsigned i = 0; i < 3; ++i) r[i] = xyz[i] < -0.5 ?
		xyz[i] + 1.0 : (xyz[i] > 0.5 ? xyz[i] - 1.0 : xyz[i]);
}

static float vec_dot (float const q[3], float const r[3]) {
	float ret = 0.0;
	for (unsigned i = 0; i < 3; ++i) ret += q[i] * r[i];
	return ret;
}

static void mat_mult (float const m[3][3], float const x[3], float ans[3]) {
	for (unsigned i = 0; i < 3; ++i) ans[i] = vec_dot (m[i], x);
}

void mat_mult_add (float const m[3][3], float const x[3], float const y[3], float ans[3]) {
	for (unsigned i = 0; i < 3; ++i) ans[i] = vec_dot (m[i], x) + y[i];
}

static float vol_get (float const angles[3]) {
	float cs[3];
	for (unsigned i = 0; i < 3; ++i) cs[i] = cos (angles[i] * PI / 180.0);
	return sqrt (1 + 2 * cs[0] * cs[1] * cs[2] - vec_dot (cs, cs));
}
// show m*n
static void d_mat_mult (float trans[3][3], float const m[3][3], float const n[3][3])	{
	float temp[3] = {0};
	for (int i=0; i<3; i++){
		for (int t=0; t<3; t++) temp[t] = 0;
		for (int j=0; j<3; j++){
			for (int k=0; k<3; k++) temp[k] += m[i][j]*n[j][k];
		}
		for (int k=0; k<3; k++){
			trans[i][k] = temp[k];
		}
	}
}

// the unit of angles is degree 
void rot_mat_mk (float trans[3][3],float const abc[3], float const angles[3], float dcm_angles[3]){
	float cs[3], vol, s,rot_cs[3], rot_sn[3];
	unsigned i;
	for (i = 0; i < 3; ++i){
		cs[i] = cos (angles[i] * PI / 180.0);
		rot_cs[i] = cos (dcm_angles[i] * PI / 180.0);
		rot_sn[i] = sin (dcm_angles[i] * PI / 180.0);
	}
	vol = vol_get(angles);
	s = sin (angles[2] * PI / 180.0);
	float car_to_fra[3][3] = {	
		{1/abc[0], -cs[2] / (abc[0]*s), (cs[0]*cs[2]-cs[1]) / (abc[0]*vol*s) },
		{0.0,    1/(abc[1] * s),     (cs[1]*cs[2]-cs[0]) / (abc[1]*vol*s) },
		{0.0,    0.0,	s/(abc[2] * vol) }
	},
	dcm_rot[3][3] = {
		{rot_cs[0]*rot_cs[1], rot_cs[0]*rot_sn[1]*rot_sn[2]-rot_sn[0]*rot_cs[2], rot_cs[0]*rot_sn[1]*rot_cs[2]+rot_sn[0]*rot_sn[2]},
		{rot_sn[0]*rot_cs[1], rot_sn[0]*rot_sn[1]*rot_sn[2]+rot_cs[0]*rot_cs[2], rot_sn[0]*rot_sn[1]*rot_cs[2]-rot_cs[0]*rot_sn[2]},
		{-rot_sn[1], rot_cs[1]*rot_sn[2], rot_cs[1]*rot_cs[2]}
	};
	d_mat_mult(trans, car_to_fra, dcm_rot);
}

struct metric *metric_mk (float const abc[3], float const angles[3]) {
	struct metric *met;
	float cs[3], r[3], s;
	unsigned i, j;
	if (!(met = malloc (sizeof (struct metric)))) return NULL;

	for (i = 0; i < 3; ++i) cs[i] = cos (angles[i] * PI / 180.0);
	s = sin (angles[2] * PI / 180.0);
	memcpy (met->trans, (float []) {
		abc[0], abc[1] * cs[2], abc[2] * cs[1],
		0.0,    abc[1] * s,     abc[2] * (cs[0] - cs[1] * cs[2]) / s,
		0.0,    0.0,            abc[2] * vol_get (angles) / s
	}, 9 * sizeof (float));
	for (i = 0; i < 8; ++i) {
		for (j = 0; j < 3; ++j) r[j] = (i >> j) & 0x1;
		mat_mult (met->trans, r, met->vert8[i]);
	}

#ifdef BENCHMARK
	float q[3];
	unsigned k;
	for (i = 0; i < 27; ++i) {
		for (j = 0, k = i; j < 3; ++j, k /= 3) r[j] = (float) (k % 3) - 1.0;
		mat_mult (met->trans, r, met->vert27[i]);
	}
	memcpy (met->g, (float []) {
		abc[0] * abc[0], abc[0] * abc[1] * cs[2], abc[0] * abc[2] * cs[1],
		abc[0] * abc[1] * cs[2], abc[1] * abc[1], abc[1] * abc[2] * cs[0],
		abc[0] * abc[2] * cs[1], abc[1] * abc[2] * cs[0], abc[2] * abc[2]
	}, 9 * sizeof (float));
	met->sq[0] = 0.0;
	for (i = 1; i < 8; ++i) {
		for (j = 0; j < 3; ++j) r[j] = (i >> j) & 0x1;
		mat_mult (met->g, r, q);
		met->sq[i] = vec_dot (q, r);
	}
#endif
	return met;
}

// <https://github.com/krishkshir/crystalLattice/blob/master/CVP.cpp>.
float metric_get_fast (struct metric const *met, float const xyz[3]) {
	float r0[3], r1[3];
	vec_round (xyz, r0);
	mat_mult (met->trans, r0, r1);
	return vec_dot (r1, r1);
}

float metric_get_slow (struct metric const *met, float const xyz[3]) {
	float r0[3], r1[3], sq[8];
	unsigned i, j;
	for (i = 0; i < 3; ++i) r1[i] = xyz[i] < 0.0 ? xyz[i] + 1.0 : xyz[i];
	mat_mult (met->trans, r1, r0);
	for (i = 0; i < 8; ++i) {
		sq[i] = 0.0;
		for (j = 0; j < 3; ++j) {
			r1[j] = r0[j] - met->vert8[i][j];
			sq[i] += r1[j] * r1[j];
		}
	}
	for (i = 1, j = 0; i < 8; ++i) if (sq[i] < sq[j]) j = i;
	return sq[j];
}

float metric_get (struct metric const *met, float const xyz[3]) {
	float r0[3], r1[3];
	vec_round (xyz, r0);
	if (fabs (r0[0]) + fabs (r0[1]) + fabs (r0[2]) < 0.5) {
		mat_mult (met->trans, r0, r1);
		return vec_dot (r1, r1);
	} else return metric_get_slow (met, xyz);
}

#ifdef BENCHMARK
float metric_get_fast_bm (struct metric const *met, float const xyz[3]) {
	float q[3], r[3];
	vec_round (xyz, r);
	mat_mult (met->g, r, q);
	return vec_dot (q, r);
}

float metric_get_slow_bm (struct metric const *met, float const xyz[3]) {
	float q[3], r[3], sq[8];
	unsigned i, j;

	for (i = 0; i < 3; ++i) r[i] = xyz[i] < 0.0 ? xyz[i] + 1.0 : xyz[i];
	mat_mult (met->g, r, q);
	sq[0] = vec_dot (q, r);

	for (i = 1; i < 8; ++i) {
		sq[i] = sq[0] + met->sq[i];
		for (j = 0; j < 3; ++j) if ((i >> j) & 0x1) sq[i] -= 2.0 * q[j];
	}
	for (i = 1, j = 0; i < 8; ++i) if (sq[i] < sq[j]) j = i;
	return sq[j];
}

float metric_get_brute (struct metric const *met, float const xyz[3]) {
	float r[3], sq[27];
	unsigned i, j;
	mat_mult (met->trans, xyz, r);
	for (i = 0; i < 27; ++i) {
		sq[i] = 0.0;
		for (j = 0; j < 3; ++j) {
			float x = r[j] + met->vert27[i][j];
			sq[i] += x * x;
		}
	}
	for (i = 1, j = 0; i < 27; ++i) if (sq[i] < sq[j]) j = i;
	return sq[j];
}

void ratios_get (float const abc[3], float const angles[3], float ratios[3]) {
	float v = vol_get (angles);
	for (unsigned i = 0; i < 3; ++i) ratios[i] =
		sin (angles[i] * PI / 180.0) / (v * abc[i]);
}
#endif

