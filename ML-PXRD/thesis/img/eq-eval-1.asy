import utils;
size (8cm);

pen c[] = {red, bleach (blue, 0.66), grey};
pair s = (2.5, 2), a = (0, 0.75), b = (1.5, 1.1);
real r[] = {0.5, 0.6};

void atom (pair pos, real rad, string lab, int i) {
	filldraw (circle (pos, rad), bleach (c[i], 0.66), c[i]);
	label ("\large" + lab, pos);
}

draw (s -- (s.x, -s.y) -- -s -- (-s.x, s.y) -- cycle, c[2]);
draw ((-s.x, 0) -- (s.x, 0), dashed + c[2]);
draw ((0, -s.y) -- (0, s.y), dashed + c[2]);

atom (a, r[0], "00", 0);
atom (-a, r[0], "01", 1);

atom (b, r[1], "10", 0);
atom ((-b.x, b.y), r[1], "11", 1);
atom (-b, r[1], "12", 1);
atom ((b.x, -b.y), r[1], "13", 1);

