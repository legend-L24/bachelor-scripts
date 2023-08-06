#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
#include "cryst.h"
#include "optim.h"
#include "utils.h"

#define ERR_GOTO(DESC, DEST) { desc = DESC; goto DEST; }
int main (int argc, char const *const *argv) {
	rng_t *rng;
        crystal *cr;
        char const *desc = NULL;
        int ret = 1;
        if (!(rng = rng_mk2 ())) ERR_GOTO ("failed rng_mk2()", rng_err)
        if (!(cr = cryst_read (rng, stdin))) ERR_GOTO ("failed cryst_read()", cr_err)
        if (!init_scan_var (cr, argc, argv)) ERR_GOTO ("failed init_scan_var()", retn_err)
        if (!mc_bias(cr)) ERR_GOTO ("failed mc_bay", retn_err)
        ret = 0; retn_err:
        cryst_fin (cr); cr_err:
        free (rng); rng_err:
        if (desc) fprintf (stderr, "E: %s\n", desc);
        return ret;
}

