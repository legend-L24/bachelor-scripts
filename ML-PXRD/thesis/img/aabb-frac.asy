import utils;
picture p[] = unit_pics (2, 3cm);
pen c[] = {red, bleach (blue, 0.66), grey};
pair xy[] = {(0.4, 0.4), (0.7, 0.7)};
real abt[] = {1.3, 1, 3 / 8 * pi}, r[] = {0.3, 0.2};

transform t[] = {(0, 0, abt[0], abt[1] * cos(abt[2]), 0, abt[1] * sin(abt[2]))};
t.push (inverse (t[0]));
t.push (scale (1 / (abt[0] * sin(abt[2])), 1 / (abt[1] * sin(abt[2]))));

draw (p[0], t[0] * unitsquare, c[2]);
draw (p[1], unitsquare, c[2]);
for (int i = 0; i < 2; ++i) {
	path atom = circle (t[0] * xy[i], r[i]);
	filldraw (p[0], atom, bleach(c[1], 0.66), c[1]);
	filldraw (p[1], t[1] * atom, bleach(c[1], 0.66), c[1]);
}
for (int i = 0; i < 2; ++i) {
	path bb = shift (xy[i]) * scale (r[i]) * t[2] * usq;
	draw (p[0], t[0] * bb, c[0]);
	draw (p[1], bb, c[0]);
}

attach (p[0].fit (), (0, 0), W);
attach (p[1].fit (), (0.6, 0), E);
draw ((0.1, 0) -- (0.4, 0), EndArrow (6));

