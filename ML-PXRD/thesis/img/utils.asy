path usq = scale(2) * shift (-0.5, -0.5) * unitsquare;
pen bleach (pen p, real x) { return x * p + (1 - x) * white; }

picture[] unit_pics (int n, real s) {
	picture p[];
	for (int i = 0; i < n; ++i) {
		p.push (new picture);
		unitsize (p[i], s);
	}
	unitsize (s);
	return p;
}

real[][] abr (real x[][]) {
	real r[][];
	for (int i = 0; i < x.length; ++i) r.push (new real[] {
		(x[i][1] + x[i][0]) / 2, (x[i][3] + x[i][2]) / 2,
		(x[i][1] - x[i][0]) / 2, (x[i][3] - x[i][2]) / 2
	});
	return r;
}

bool edge (real b[], real x) {
	for (int i = 0; i < b.length; ++i) {
		if (abs (x - b[i]) < 0.01) return true;
	}
	return false;
}

void ticks (real b[], real x[][]) {
	for (int i = 0; i < x.length; ++i) {
		for (int j = 0; j < 2; ++j) if (!edge (b, x[i][j])) {
			draw ((x[i][j], 0) -- (x[i][j], x[i][2]), dashed);
		}
	}
}

void aabbs (pen c, real b[], real x[][]) {
	for (int i = 0; i < x.length; ++i) {
		for (int j = 0; j < 2; ++j) {
			draw ((x[i][0], x[i][2 + j]) -- (x[i][1], x[i][2 + j]), c);
			if (!edge (b, x[i][j])) draw
				((x[i][j], x[i][2]) -- (x[i][j], x[i][3]), c);
		}
	}
	for (int i = 0; i < x.length; ++i) label
		("\large" + string (x[i][6]), (x[i][4], x[i][5]));
}

string fbox (string s) {
	return "\fbox{" + s + "}";
}

string tag (string w, string s) {
	return "\parbox{" + w + "}{\centering " + s + "}";
}

