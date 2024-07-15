#include <stdlib.h>
#include "dftbplus_plugin.h"

int counter = -1;
int nAtoms = 0;

int init() {
    counter = 0;
    return 1;
}

void setNeighbourList(int n, int nAtomsCent, double *coords, int *img2CentCell, int *iNeighbour, double *neightDist2) {
    nAtoms = n;
}


int getSKIntegrals(int nSkgrid, int nSkIntg, double *skTab, double dist,
    int atom1, int atom2, int species1, int species2, int HorS, double interdist) {
    skTab[0] = 0.0;
    if (nAtoms != 2 || species1 != 1 || species2 != 1) {
        return 0;
    }

    if (counter == 0) {
        skTab[0] = 1.0;
    } else if (counter == 1) {
        skTab[0] = 2.0;
    }
    counter++;
    return 1;
}
