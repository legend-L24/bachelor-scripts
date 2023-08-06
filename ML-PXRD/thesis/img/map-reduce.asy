unitsize (2.4cm, 1.6cm);
arrowbar a = EndArrow (5);
string s[][] = {
	{"$x_0$", "$x_1$", "$\cdots$", "$x_n$"},
	{"$f(x_0)$", "$f(x_1)$", "$\cdots$", "$f(x_n)$"}
};

label ("$\bigoplus_i f(x_i)$", (0, 2));
for (int i = 0; i < 4; ++i) {
	real x = i - 1.5;
	label (s[0][i], (x, 0));
	label (s[1][i], (x, 1));
	draw ((x, 0.25) -- (x, 0.75), a);
	draw ((x, 1.25) -- (x / 8, 1.75), a);
}

