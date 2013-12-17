#include "../include/Profile.h"

/*
Created:        13 december 2013    by Jeroen Barnhoorn
Last updated:   13 december 2013    by Jeroen Barnhoorn

STATUS: DONE
*/

Profile::Profile() {
}

GridBdry Profile::getGridBdry() {
    GridBdry gridbdry(Gamma_I, interior, n_grid);
    return gridbdry;
}

ExactBdry Profile::getExactBdry() {
    ExactBdry exactbdry(ex_loc, ex_ds, ex_n);
    return exactbdry;
}
