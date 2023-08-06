extern float metric_get_fast_bm (metric const *met, float const xyz[3]);
extern float metric_get_slow_bm (metric const *met, float const xyz[3]);
extern float metric_get_brute (metric const *met, float const xyz[3]);
extern void ratios_get
	(float const abc[3], float const angles[3], float ratios[3]);

