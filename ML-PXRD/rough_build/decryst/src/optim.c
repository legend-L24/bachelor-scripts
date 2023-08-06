#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"
#include "cryst.h"
#include "optim.h"

struct sa_comp {
	crystal *cr;
	rng_t *rng;
	step_t *step;
	uint32_t *bbuf;
	float best[3];
};

struct sa_sched {
	// f: (alpha, f(1), f(s), f(s^2), f(1/u), f(s/u));
	// g: (beta, g(1), g(s), g(s^2), g(1/v), g(s/v)).
	float f[6], g[6], mag, lambda;
};

int mc_stat (crystal *cr, rng_t *rng, unsigned n, float ret[5]) {
	if (!cryst_eval (cr)) return 0;
	cryst_ack (cr, 1);
	float const *scores = cryst_scores (cr);
	memcpy (ret + 2, scores, 3 * sizeof (float));
	ret[0] = ret[2];
	ret[1] = ret[2] * ret[2];
	for (unsigned i = 1; i < n; ++i) {
		cryst_step (cr, rng_float (rng),1);
		if (!cryst_eval (cr)) return 0;
		cryst_ack (cr, 1);
		float x = scores[0];
		ret[0] += x;
		ret[1] += x * x;
		if (x < ret[2]) memcpy (ret + 2, scores, 3 * sizeof (float));
	}
	ret[0] /= n;
	ret[1] = (ret[1] - ret[0] * ret[0] / n) / (n - 1);
	if (ret[1] < 1e-12) ret[1] = 1e-12;
	ret[1] = sqrt (ret[1]);
	return 1;
}

struct sa_sched *sa_sched_mk (
	float const ab[2], float const init[3], float ctl[7]
) {
	struct sa_sched *sched;
	if (!(sched = malloc (sizeof (struct sa_sched)))) return NULL;
	*sched = (struct sa_sched) { .mag = 1.0, .lambda = init[0] };
	memcpy (sched->f, (float []) {
		ab[0], 1.0, 0.0, 0.0, 1.0 / init[1], 0.0
	}, 6 * sizeof (float));
	memcpy (sched->g, (float []) {
		ab[1], 1.0, 0.0, 0.0, 1.0 / init[2], 0.0
	}, 6 * sizeof (float));
	memcpy (ctl, (float []) {
		0.0, sched->f[4], init[1] / init[2], sched->g[4],
		sched->mag, 0.0, 0.0
	}, 7 * sizeof (float));
	ctl[0] = ctl[2] * ctl[2];
	return sched;
}

static void ab_get (float const arr[6], float ret[2]) {
	if (!isfinite (
		ret[0] = (arr[1] * arr[5] - arr[2] * arr[4]) /
		(arr[1] * arr[3] - arr[2] * arr[2])
	)) ret[0] = (arr[1] * arr[5] + arr[2] * arr[4]) /
		(arr[1] * arr[3] + arr[2] * arr[2]);
	ret[1] = (arr[4] - ret[0] * arr[2]) / arr[1];
}

int sa_sched_step (struct sa_sched *sched, float const stat[4], float ret[6]) {
	float now[6] = {
		0.0, 1.0, stat[2], stat[2] * stat[2],
		1.0 / stat[0], stat[2] / stat[0]
	};
	unsigned i;

	for (i = 1; i < 6; ++i) sched->f[i] = sched->f[0] * sched->f[i] + now[i];
	now[4] = 1.0 / stat[1];
	now[5] = stat[2] / stat[1];
	for (i = 1; i < 6; ++i) sched->g[i] = sched->g[0] * sched->g[i] + now[i];
	ab_get (sched->f, ret);
	ab_get (sched->g, ret + 2);
	for (i = 0; i < 4; ++i) if (!isfinite (ret[i])) return 0;

	sched->mag *= stat[3] + 0.56;
	if (sched->mag < 1e-6) sched->mag = 1e-6;
	else if (sched->mag > 1.0) sched->mag = 1.0;
	ret[4] = sched->mag;

	ret[5] = (1.0 - stat[3]) / (2.0 - stat[3]);
	ret[5] *= 4.0 * sched->lambda * stat[3] * ret[5];
	return 1;
}

struct sa_comp *sa_comp_mk (crystal *cr, rng_t *rng) {
	struct sa_comp *comp;
	if (!(comp = malloc (sizeof (struct sa_comp)))) goto comp_err;
	*comp = (struct sa_comp) { .cr = cr, .rng = rng, .best = { 1e9, 0.0 } };
	if (!(comp->bbuf = calloc (
		cryst_dof (cr), sizeof (uint32_t))
	)) goto bbuf_err;
	if (!(comp->step = step_mk (rng))) goto step_err;
	return comp;
	free (comp->step); step_err:
	free (comp->bbuf); bbuf_err:
	free (comp); comp_err:
	return NULL;
}

void sa_comp_fin (struct sa_comp *comp) {
	free (comp->step);
	free (comp->bbuf);
	free (comp);
}

float const *sa_comp_best (struct sa_comp *comp) {
	return comp->best;
}

uint32_t const *sa_comp_bbuf (struct sa_comp *comp) {
	return comp->bbuf;
}

int sa_comp_get (
	struct sa_comp *comp, unsigned const n[3], float const ctl[8],
    unsigned iret[1], float fret[3], unsigned flag
) {
	float const *scores = cryst_scores (comp->cr);
	float fr[2] = { 0.0, 0.0 }, x = scores[0], r = ctl[5], s = ctl[6], tmp;
	unsigned cnt = 0;

	if (s > 0.0) {
		tmp = ctl[2] * s + ctl[3];
		s += (n[1] + 1) * r / (s * s) * tmp * tmp * tmp;
	} else {
		s = 0.5 * ctl[3];
		tmp = ctl[2] * s + ctl[3];
		s += n[1] * r / (s * s) * tmp * tmp * tmp;
	}
	r *= n[2];
	step_mag (comp->step, ctl[4]);

	for (unsigned i = 0;;) {
		cryst_step (comp->cr, step_get (comp->step), flag);
		if (!cryst_eval (comp->cr)) return 0;
		float y = scores[0];
		if (y < comp->best[0]) {
			memcpy (comp->best, scores, 3 * sizeof (float));
			cryst_dump (comp->cr, comp->bbuf);
		}

		tmp = y - x;
		if (tmp > 0.0 && rng_float (comp->rng) > exp (-s * tmp)) {
			cryst_ack (comp->cr, 0);
		} else {
			cryst_ack (comp->cr, 1);
			x = y;
			++cnt;
		}

		fr[0] += x;
		tmp = x - 1.0 / (ctl[0] * s + ctl[1]);
		if (!isfinite (tmp)) tmp = x - fr[0] / (i + 1);
		fr[1] += tmp * tmp;

		if (++i < n[0]) {
			tmp = ctl[2] * s + ctl[3];
			tmp = s + r / (s * s) * tmp * tmp * tmp;
			s *= ctl[7];
			if (tmp < s) s = tmp;
		} else break;
	}

	memcpy (fret, fr, 2 * sizeof (float));
	fret[2] = s;
	iret[0] = cnt;
	return 1;
}

