import utils;
texpreamble ("\usepackage{isomath}");
picture p = unit_pics (1, 1.2cm)[0];
pen d = linetype (new real[] {4,4}), g = 0.33 * white;
arrowbar a = EndArrow (4);

string lab (string cnt, int i) {
	return "$\tensorsym{x}_{" + (
		i == 2 ? "\cdots" : string (i)
	) + "," + cnt + "}$";
}

for (int i = 0; i < 3; ++i) {
	int j = 1 - i;
	draw (p, (-3.5, 0.6 * j) -- (-3, 0.6 * j), a);
	draw (p, (-1, 0.6 * j) -- (-0.5, 0.6 * j));
	draw (p, (-0.5, 0.6 * j) -- (0.5, 0.6 * j), dotted);
	draw (p, (0.5, 0.6 * j) -- (1, 0.6 * j), a);
	draw (p, (3, 0.6 * j) -- (3.5, 0.6 * j));
	label (fbox (tag ("4.6em", lab ("n\tau", i))), (-6, 0.6 * j));
	label (fbox (tag ("4.6em", lab ("(n + 1)\tau - 1", i))), (-2, 0.6 * j));
	label (fbox (tag ("4.6em", lab ("(n + 1)\tau", i))), (2, 0.6 * j));
	label (fbox (tag ("4.6em", lab ("(n + 2)\tau - 1", i))), (6, 0.6 * j));
}

draw (p, (-3.5, -0.6) -- (-3.5, 0.6));
draw (p, (3.5, -0.6) -- (3.5, 0.6));
draw (p, (-0.2, 0.5) -- (-0.2, -1), d + g, a);
draw (p, (0.2, -0.7) -- (0.2, -1), d + g, a);
draw (p, (0, -0.1) -- (0, -1), d + g, a);
label (p, "moves made in parallel", (0, -1), S, g);

draw ((-8, 0) -- (-7.5, 0));
draw ((-0.5, 0) -- (0.5, 0));
draw ((7.5, 0) -- (8, 0), a);
draw ((0, 0.1) -- (0, 1), d + g, a);
label ("update of estimators and (every several periods) mixing of states", (0, 1), N, g);
label ("$\cdots$", -8, W);
label ("$\cdots$", 8, E);
attach (p.fit (), (-4, 0));
attach (p.fit (), (4, 0));

