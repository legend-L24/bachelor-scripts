import graph;
size (9cm, 6cm, keepAspect = false);

draw ((0, 0) -- (2, 0), grey, EndArrow (size = 6));
draw ((0, -1.25) -- (0, 1.25), grey, EndArrow (size = 6));
draw ((1, 0) -- (1, -1), dashed + grey);
draw (graph (new real (real x) {
	return 1 / x ^ 12 - 2 / x ^ 6;
}, 0.855, 1.8));
draw ((0, 1) -- (0.75, 1) -- (0.875, 0) -- (1.8, 0), red);

label ("0", (0, 0), W);
label ("$d / d_0$", (2, 0), E);
label ("$f (d / d_0)$", (0, 1.25), N);
label ("1", (1, 0), N);

