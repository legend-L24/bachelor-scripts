unitsize (1.2cm);
texpreamble ("\usepackage{xcolor}");
texpreamble ("\newcommand{\ul}{\string_}");
texpreamble (
	"\newcommand{\head}[2]{\parbox{10em}{\strut#1\\"
	"\textcolor{black!60}{\footnotesize\textit{\strut#2}}}}"
);

arrowbar a = EndArrow (5);
string desc (string s) { return "\parbox{7em}{\footnotesize{" + s + "}}"; }
void part (real lbt[], string head[], string body[]) {
	real h[] = {lbt[2] - 0.6, lbt[2]},
		lrbt[] = {lbt[0] - 0.1, lbt[0] + 3.5, lbt[1], lbt[2] + 1};
	label ("\head{" + head[0] + "}{" + head[1] + "}", (lbt[0], lbt[2]), NE);
	for (int i = 0; i < body.length; ++i) label (
		"\texttt{{\strut}" + replace (body[i], "_", "\string_") + "}",
		(lbt[0], h[0] - 0.5 * i), NE
	);
	draw (
		(lrbt[0], lrbt[2]) -- (lrbt[1], lrbt[2]) --
		(lrbt[1], lrbt[3]) -- (lrbt[0], lrbt[3]) -- cycle
	);
	draw ((lrbt[0], h[1]) -- (lrbt[1], h[1]));
}

part (
	new real[] {0, 5.3, 7},
	new string[] {"Control program", "main user interface"},
	new string[] {"decr_utils"}
);
part (
	new real[] {0, 0, 2.2},
	new string[] {"Helper scripts", "glues and abstractions"},
	new string[] {
		"do_stat.sh", "do_filter.sh",
		"do_optim.sh", "\strut$\cdots$"
	}
);
part (
	new real[] {6.6, 5.3, 7},
	new string[] {"Wrapper scripts", "distributed execution"},
	new string[] {"stat.sh", "optim.sh", "server.sh"}
);
part (
	new real[] {6.6, 0, 2.2},
	new string[] {"Lower-level programs", "fundamental tasks"},
	new string[] {"decr_mcs", "decr_lsa", "decr_sas", "decr_sac"}
);

label (desc ("abstract and\\bridge some\\common tasks"), (2.1, 4.25));
label (desc ("runs wrappers\\in parallel\\as multiple\\processes"), (5.3, 6.4));
label (desc ("run lower-level\\programs on\\specified machines"), (8.7, 4.25));
draw ((0.5, 3.25) -- (0.5, 5.25), a);
draw ((3.55, 7.4) -- (6.45, 7.4), a);
draw ((7.1, 5.25) -- (7.1, 3.25), a);

