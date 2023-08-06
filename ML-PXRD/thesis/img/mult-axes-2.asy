size (8cm);
draw (scale (3, 5) * unitsquare, grey);

for (int i = 0; i < 3; ++i) draw (circle ((0.9, i + 1.5), 0.6));
for (int i = 0; i < 4; ++i) draw (circle ((2.1, i + 1), 0.6));

draw ((0.25, 0.15) -- (0.25, 0.5), EndArrow (5));
draw ((0.25, 0.15) -- (0.6, 0.15), EndArrow (5));
label ("$x$", (0.6, 0.15), E);
label ("$y$", (0.25, 0.5), N);

