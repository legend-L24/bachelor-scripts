import utils;
size (0, 4cm);
pen c[] = {red, bleach (blue, 0.5), grey};
real h = 1.9, r = 0.025, g = 1 / (8 * h);

fill ((1, 2 * h / 3) -- (1, h / 2 + g) -- (0.75, h / 2) -- cycle, c[1]);
fill ((1.5, h / 3) -- (1.5, h / 2 - g) -- (1.75, h / 2) -- cycle, c[1]);
for (int i = 0; i < 3; ++i) {
	draw (
		(i - 0.25, h / 2) -- (i, h / 2 + g) --
		(i + 0.5, h / 2 - g) -- (i + 0.75, h / 2), c[2]
	);
	draw ((i, h / 2 + g) -- (i, h + 2 * g), c[2]);
	draw ((i + 0.5, h / 2 - g) -- (i + 0.5, - 2 * g), c[2]);
}
draw ((0, 0) -- (1, 0) -- (2.5, h) -- (1.5, h) -- cycle);
for (int i = 0; i < 3; ++i) {
	fill (circle ((i, 0), r), c[0]);
	fill (circle ((i + 0.5, h), r), c[0]);
}

