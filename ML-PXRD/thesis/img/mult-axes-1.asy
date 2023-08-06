size (8cm);
draw (scale (3) * unitsquare, grey);

for (int i = 0; i < 2; ++i) {
	for (int j = 0; j < 2; ++j) draw (circle ((i + 1, j + 1), 0.4));
	for (int j = 0; j < 3; ++j) {
		draw (circle ((i + 1, j + 0.5), 0.2));
		draw (circle ((j + 0.5, i + 1), 0.2));
	}
}

draw ((0.15, 0.15) -- (0.15, 0.4), EndArrow (5));
draw ((0.15, 0.15) -- (0.4, 0.15), EndArrow (5));
label ("$x$", (0.4, 0.15), E);
label ("$y$", (0.15, 0.4), N);

