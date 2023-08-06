import utils;
unitsize (1.2cm);
texpreamble ("\newcommand*{\mytimes}{\!\times\!}");
arrowbar a = EndArrow (4);
pen g = 0.33 * white;
string t[] = {
	"Pb$^{2+}$: $4a\mytimes1$; S$^{6+}$: $4b\mytimes1$; O$^{2-}$: $4c\mytimes4$",
	"Pb$^{2+}$: $4a\mytimes1$; S$^{6+}$: $4b\mytimes1$; O$^{2-}$: $4c\mytimes2\!+\!8d\mytimes1$",
	"Pb$^{2+}$: $\cdots$; S$^{6+}$: $\cdots$; O$^{2-}$: $\cdots$"
};

label (fbox (tag ("4em", "PbSO$_4$ ($Pnma$)")), (0.05, 0), W);
label (tag ("5em", "tasks done in parallel"), (9.9, 0), E, g);
draw (brace ((9.7, 0.9), (9.7, -0.9)), g);
for (int i = 0; i < 3; ++i) {
	int j = 1 - i;
	label (fbox (tag ("19.5em", t[i])), (4.45, 0.8 * j));
	label ("$\cdots$", (9.3, 0.8 * j));
	draw ((0.1, 0.4 * j) -- (0.7, 0.8 * j), a);
	draw ((8.2, 0.8 * j) -- (8.9, 0.8 * j), a);
}

