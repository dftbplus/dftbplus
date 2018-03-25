#ifndef __DFTBPLUS_H__
#define __DFTBPLUS_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  void *pDftbPlusInput;
} DftbPlusInput;
  
typedef struct DftbPlus {
  void * pDftbPlus;
} DftbPlus;

void dftbp_init(DftbPlus *instance);
void dftbp_destruct(DftbPlus *instance);
void dftbp_get_input_from_file(DftbPlus *instance, const char *filename, DftbPlusInput *input);
void dftbp_setup_calculator(DftbPlus *instance, DftbPlusInput *input);
void dftbp_set_coords(DftbPlus *instance, double *coords);
void dftbp_set_coords_and_latvecs(DftbPlus *instance, double *coords, double *latvecs);
void dftbp_get_energy(DftbPlus *instance, double *mermin_energy);
void dftbp_get_gradients(DftbPlus *instance, double *gradients);
  
#ifdef __cplusplus
}
#endif

#endif
