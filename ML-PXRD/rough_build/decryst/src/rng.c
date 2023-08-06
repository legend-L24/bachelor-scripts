// Minimal version of the PCG RNG <http://pcg-random.org/>, chosen for its
// simplicity, high speed, not quite bad statistical characteristics and not too
// short period (2^64, which the author of this software considers to be quite
// adequate for most crystallographical applications).

/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *       http://www.pcg-random.org
 */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "rng.h"

#define SEED_FILE "/dev/urandom"

struct rng_t {
	uint64_t state, inc;
};

struct step_t {
	struct rng_t *rng;
	float a, b;
};

struct rng_t *rng_mk (uint64_t seed) {
	struct rng_t *rng;
	if (!(rng = malloc (sizeof (struct rng_t)))) return NULL;
	*rng = (struct rng_t) { .state = 0, .inc = seed << 1 | 0x1 };
	rng_next (rng);
	rng->state += seed;
	rng_next (rng);
	return rng;
}

uint32_t rng_next (struct rng_t *rng) {
	uint32_t rot = rng->state >> 59,
		xorshifted = ((rng->state >> 18) ^ rng->state) >> 27;
	rng->state = rng->state * UINT64_C (0x5851f42d4c957f2d) + rng->inc;
	return xorshifted >> rot | xorshifted << ((-rot) & 0x1f);
}

double rng_float (struct rng_t *rng) {
	return rng_next (rng) * 2.32830643653869628906e-10;
}

unsigned rng_dice (struct rng_t *rng, unsigned n) {
	uint32_t drop = -n;
	while (1) {
		uint32_t x = rng_next (rng), ret = x % n;
		if (x - ret <= drop) return ret;
	}
}

unsigned rng_bdice (struct rng_t *rng, float const qs[], unsigned n) {
	float x = rng_float (rng);
	unsigned a = 0, b = n, i;
	while (a < b) {
		i = a + (b - a - 1) / 2;
		if (x == qs[i]) return i;
		else if (x < qs[i]) b = i;
		else a = i + 1;
	}
	return a;
}

void rng_shuf (struct rng_t *rng, unsigned arr[], unsigned n) {
	unsigned i, j, tmp;
	for (i = n - 1; i; --i) {
		j = rng_dice (rng, i + 1);
		tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
	}
}

// Wrapped exponential distribution, scaled and with random sign.
struct step_t *step_mk (struct rng_t *rng) {
	struct step_t *step;
	if (!(step = malloc (sizeof (struct step_t)))) return NULL;
	step->rng = rng;
	return step;
}

void step_mag (struct step_t *step, float mag) {
	step->a = 0.5 * mag;
	step->b = exp (-1.0 / mag) - 1.0;
}

float step_get (struct step_t const *step) {
	// `x' used to enforce evaluation order of rng_*().
	float x = (rng_next (step->rng) % 2 ? 1 : -1) * step->a;
	return x * log1p (step->b * rng_float (step->rng));
}

int seed_mk (uint64_t *seed) {
	int fd;
	if ((fd = open (SEED_FILE, O_RDONLY)) == -1) goto file_err;
	if (read (fd, seed, sizeof (uint64_t)) == -1) goto read_err;
	if (close (fd) == -1) goto file_err;
	return 1; read_err:
	close (fd); file_err:
	return 0;
}

struct rng_t *rng_mk2 (void) {
	uint64_t seed;
	return seed_mk (&seed) ? rng_mk (seed) : NULL;
}

void wstoqs (float qs[], unsigned n) {
	float tmp = 0.0;
	unsigned i;
	for (i = 0; i < n; ++i) tmp += qs[i];
	for (i = 0; i < n; ++i) qs[i] /= tmp;
	for (i = 1; i < n; ++i) qs[i] += qs[i - 1];
	// Handle error accumulation in float addition.
	if (qs[n - 1] < 1.0) qs[n - 1] = 1.0;
}

