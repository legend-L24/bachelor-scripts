import utils;
size (8cm);
transform t = (0, 0, 1, 0.5, 0, sqrt (3) / 2);

for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
	draw (circle (t * (i, j), 0.5), grey);
}
draw (t * unitsquare);
draw (t * shift (1, 0) * scale (1 / sqrt (3)) * usq, red);

