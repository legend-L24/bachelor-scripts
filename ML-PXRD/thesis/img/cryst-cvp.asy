import utils;
picture p[] = unit_pics (2, 2cm);
arrowbar a = EndArrow (5);
pen c[] = {red, bleach (blue, 0.66), grey};
real d[] = {0.8, 0.2}, r[] = {0.1, 0.15}, e[] = {0.05, 0.2};
d.push (d[0] - d[1]);

void atom (picture pic, pair move, int i) {
	filldraw (pic, circle (move, r[i]), bleach (c[i], 0.66), c[i]);
}

for (int i = 0; i <= 2; ++i) {
	draw (p[0], (-e[0], i) -- (2 + e[0], i), c[2]);
	draw (p[0], (i, -e[0]) -- (i, 2 + e[0]), c[2]);
}
for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
	for (int k = 0; k < 2; ++k) atom (p[0], (i + d[k], j + d[k]), k);
	draw (p[0], (d[0], d[0]) -- (i + d[1], j + d[1]), a);
}

atom (p[1], (d[2], d[2]), 0);
for (int i = 0; i < 2; ++i) {
	draw (p[1], (-e[1], i) -- (1 + e[1], i), c[2]);
	draw (p[1], (i, -e[1]) -- (i, 1 + e[1]), c[2]);
}
for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
	atom (p[1], (i, j), 1);
	draw (p[1], (d[2], d[2]) -- (i, j), a);
}

attach (p[0].fit (), (0, 0), W);
attach (p[1].fit (), (0.8, 0), E);
draw ((0.3, 0) -- (0.7, 0), a);

