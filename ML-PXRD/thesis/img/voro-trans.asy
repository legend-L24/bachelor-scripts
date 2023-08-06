import utils;
picture p[] = unit_pics (4, 3cm);
real x[] = {0.08, 0.4}, r = 0.5;
x.push (x[0] * (0.25 - x[0]) / x[1]);
x.push (x[0] * x[1] / (0.25 + x[0]));

void dot_cell (picture pic, path pt) {
	draw (pic, pt);
	dot (pic, pt, red);
}

fill (p[2], (-0.25, 0) -- (-0.25 + x[0], x[2])
-- (-0.25 + x[0], x[3]) -- cycle, bleach (blue, 0.5));
fill (p[2], (0.25, 0) -- (0.25 - x[0], -x[2])
-- (0.25 - x[0], -x[3]) -- cycle, bleach (blue, 0.5));

draw (
	p[0], (-0.5, 0) -- (-0.25 - x[0], x[2]) -- (-0.25 + x[0], -x[2])
	-- (0.25 - x[0], x[2]) -- (0.25 + x[0], -x[2]) -- (0.5, 0), grey
);
draw (p[0], (-0.25 - x[0], x[2]) -- (-0.25 - x[0], r), grey);
draw (p[0], (-0.25 + x[0], -x[2]) -- (-0.25 + x[0], -r), grey);
draw (p[0], (0.25 - x[0], x[2]) -- (0.25 - x[0], r), grey);
draw (p[0], (0.25 + x[0], -x[2]) -- (0.25 + x[0], -r), grey);

draw (p[3], (-0.5, 0) -- (0.5, 0), grey);
draw (p[3], (-0.25, -r) -- (-0.25, r), grey);
draw (p[3], (0.25, -r) -- (0.25, r), grey);

attach (p[1], p[0].fit (), (0, 0));
attach (p[2], (xscale (-1) * p[0]).fit (), (0, 0));
dot_cell (
	p[1], (x[0], -x[1]) -- (0.5 - x[0], x[1])
	-- (-x[0], x[1]) -- (-0.5 + x[0], -x[1]) -- cycle
);
dot_cell (
	p[2], (-x[0], -x[1]) -- (0.5 + x[0], x[1])
	-- (x[0], x[1]) -- (-0.5 - x[0], -x[1]) -- cycle
);
dot_cell (
	p[3], (0, -x[1]) -- (0.5, x[1])
	-- (0, x[1]) -- (-0.5, -x[1]) -- cycle
);

attach (p[1].fit (), (-1.6, 0));
attach (p[2].fit (), (1.6, 0));
attach (p[3].fit (), (0, 0));
draw ((-0.9, 0) -- (-0.7, 0), Arrows (4));
draw ((0.7, 0) -- (0.9, 0), Arrows (4));

