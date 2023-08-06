import graph3;
import utils;
picture p[] = unit_pics (3, 2cm);
real r[] = {0.6, 0.2, 0.6}, a = 6;
currentprojection = perspective (1, 1, 1);
nmesh = 17;

fill (p[0], scale (r[2] / 0.5) * shift(-0.5, -0.5) * unitsquare, bleach (blue, 0.33));
draw (p[0], (-r[2], -r[2]) -- (-r[2], r[2]), red, EndArrow (a));
draw (p[0], (r[2], -r[2]) -- (r[2], r[2]), red, EndArrow (a));
draw (p[0], (-r[2], -r[2]) -- (r[2], -r[2]), EndArrow (a));
draw (p[0], (-r[2], r[2]) -- (r[2], r[2]), EndArrow (a));

draw (p[1], circle ((0, 0), 0.25), blue);
draw (p[1], circle ((0.9, 0), 0.25), blue);
label (p[1], "\large$\times$", (0.45, 0));

draw (p[2], surface (new triple (pair t) {
	return (
		(r[0] + r[1] * cos (t.y)) * cos (t.x),
		(r[0] + r[1] * cos (t.y)) * sin (t.x),
		r[1] * sin (t.y)
	);
}, (0, 0), (2pi, 2pi), Spline), bleach (blue, 0.33));

attach (p[0].fit (), (-0.1, 0));
attach (p[1].fit (), (1.725, 0));
attach (p[2].fit (), (4.7, 0.125));
draw ((0.75, 0) -- (1.25, 0), EndArrow (a));
draw ((3.1, 0) -- (3.6, 0), EndArrow (a));

