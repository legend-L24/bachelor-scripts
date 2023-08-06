unitsize (1.2cm);
arrowbar a = EndArrow (8);
real dist[] = {5, 1.5, 1, 0.5},
	m[] = {0.5 * dist[0], 1.5 * dist[0] - dist[2] + dist[3]};

void file (pair pos, string lab) {
	draw (shift (pos) * scale (2, 0.75) * shift (-0.5, -0.5) * unitsquare);
	label ("\texttt{" + lab + "}", pos);
}

for (int i = 0; i < 3; ++i) {
	real y = (1 - i) * dist[1];
	file ((0, y), "part" + string (i) + ".c");
	file ((dist[0], y), "part" + string (i) + ".o");
	draw ("compile", (m[0] - dist[2], y) -- (m[0] + dist[2], y), N, a);
	draw ((m[1] - dist[3], y) -- (m[1], y));
}

file ((2 * dist[0], 0), "program");
draw ((m[1], -dist[1]) -- (m[1], dist[1]));
draw ("link", (m[1] - 0.3, 0) -- (1.5 * dist[0] + dist[2], 0), N, a);

