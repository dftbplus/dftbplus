#include <stdlib.h>
#include "dftbplus_plugin.h"

int counter = -1;
int nAtoms = 0;

int init() {
    counter = 0;
    return 1;
}

void setNeighbourList(int n, double *coords, int *img2CentCell) {
    nAtoms = n;
}

int getSKIntegrals(int nSk, double *sk, double dist, int atom1, int atom2, int sp1, int sp2) {
    sk[0] = 0.0;
    if (nAtoms != 2 || sp1 != 1 || sp2 != 1) {
        return 0;
    }

    if (counter == 0) {
        sk[0] = -0.08;
    } else if (counter == 1) {
        sk[0] = 0.2;
    }
    counter++;
    return 1;
}
