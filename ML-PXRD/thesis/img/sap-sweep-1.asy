import utils;
texpreamble ("\usepackage{amssymb}");
unitsize (2cm, 0.66cm);

pen c[] = {red, bleach (blue, 0.66), 0.25 * white};
real x[][] = {
	{1, 3, 1, 4, 1.5, 2.5, 0}, {2, 5, 3, 6, 3.5, 5.5, 1},
	{4, 6, 2, 5, 5.5, 3.5, 2}
}, r[][] = abr (x), b[] = {};
string t[] = {
	"$\varnothing\leftarrow0$", "$\{0\}\leftarrow1$", "$\{1\}\rightarrow0$",
	"$\{1\}\leftarrow2$", "$\{2\}\rightarrow1$", "$\varnothing\rightarrow2$"
};

for (int i = 0; i < 3; ++i) fill
	(ellipse ((r[i][0], r[i][1]), r[i][2], r[i][3]), bleach (c[1], 0.66));
for (int i = 0; i < 3; ++i) draw
	(ellipse ((r[i][0], r[i][1]), r[i][2], r[i][3]), c[1]);
ticks (b, x);
aabbs (c[0], b, x);

draw ((0.5, 0) -- (6.5, 0), EndArrow);
for (int i = 0; i < 6; ++i) label ("\small" + t[i], (i + 1, -0.75));
label (
	"\footnotesize($S$ [$\leftarrow$/$\rightarrow$] $i$:" +
	" $i$ [added to$\,/\,$deleted from] $S$)",
	(6.33, -1.7), W
);

