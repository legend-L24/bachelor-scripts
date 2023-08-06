import utils;
texpreamble ("\usepackage{amssymb}");
unitsize (2cm, 0.66cm);

pen c[] = {red, bleach (blue, 0.66), 0.5 * white};
real x[][] = {
	{0.5, 2, 3, 4, 1.25, 3.5, 0}, {3, 6.5, 3, 4, 4.75, 3.5, 0},
	{1, 4, 1, 2, 2.5, 1.5, 1}, {5, 6, 1, 2, 5.5, 1.5, 2}
}, r[][] = abr (x), b[] = {0.5, 6.5};
string t[] = {
	"$\varnothing\leftarrow1$", "$\{1\}\not\rightarrow0$",
	"$\{1\}\leftarrow0$", "$\{0\}\rightarrow1$",
	"$\{0\}\leftarrow2$", "$\{0\}\rightarrow2$",
	"$\{0\}\not\leftarrow1$", "$\varnothing\rightarrow0$"
};

ticks (b, x);
for (int i = 0; i < 4; ++i) fill (
	shift (r[i][0], r[i][1]) * scale(r[i][2], r[i][3]) * usq,
	bleach (c[1], 0.66)
);
aabbs (c[0], b, x);

draw ((0.5, 0) -- (6.5, 0));
draw ((0.5, 0) -- (0.5, 5), c[2]);
draw ((6.5, 0) -- (6.5, 5), c[2]);
for (int i = 0; i < 6; ++i) label ("\small" + t[i], (i + 1, -0.75));
for (int i = 0; i < 2; ++i) label ("\small" + t[i + 6], (i + 1, -1.7));
label (
	"\footnotesize($S$ [$\not\leftarrow$/$\not\rightarrow$] $i$:" +
	" $i$ not [added to$\,/\,$deleted from] $S$)",
	(6.4, -1.7), W
);

