size (6cm);
real t = 0.3 pi, r[][] = {{-1.5, 2}, {-1.5, 1.75}};

draw (unitcircle);
draw ((r[0][0], 0) -- (r[0][1], 0));
draw ((r[0][0], r[1][1]) -- (r[0][1], r[1][1]));
draw ((r[0][0], -1) -- (r[0][1], -1), dashed);
draw ((r[0][0], 1) -- (r[0][1], 1), dashed);
draw ((0, r[1][0]) -- (0, r[1][1]), EndArrow);
draw (
	(r[1][0] / tan (t), r[1][0]) --
	(r[1][1] / tan (t), r[1][1]), EndArrow
);

label ("$O$", (0, 0), NW);
label ("$bOc$", (r[0][1], 0), NW);
label ("$r$", (0, 1), NW);
label ("$R$", (1 / tan(t), 1), NW);
label ("$a^*$", (0, r[1][1]), N);
label ("$a$", (r[1][1] / tan(t), r[1][1]), N);

