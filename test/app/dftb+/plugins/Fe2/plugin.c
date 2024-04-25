#include <stdlib.h>
#include <stdio.h>
#include "dftbplus_plugin.h"

int counter = -1;

void getSKIntegrals(double *sk, double dist, int sp1, int sp2) {
    sk[0] = 0.0;
    if (sp1 == 1 && sp2 == 1) {
        if (counter == 0) {
            sk[0] = -0.08;
        } else if (counter == 1) {
            sk[0] = 0.2;
        }
    }
    counter++;
}

void init() {
    counter = 0;
}
