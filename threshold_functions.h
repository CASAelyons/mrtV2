/* $Id: threshold_functions.h,v 0.9 2014-12-10 11:32:23 elyons Exp $ */
/* Copyright 2014 University of Massachusetts Amherst all rights reserved*/

/* threshold_functions.h declares common functions and structs used in the UMass suite of thresholding algorithms */

struct index_pair {
  int xind;
  int yind;
};

struct latLonBrng {
  double lat;
  double lon;
  double brng;
};

struct latLon {
  double lat;
  double lon;
};

short ip_is_equal(struct index_pair a, struct index_pair b);
double degToRad(double degs);
double radToDeg(double rads);

int vincenty(double inLat, double inLon, double inDx, double inBrng, struct latLonBrng *llb_out);
int inverse_vincenty(double inLat1, double inLon1, double inLat2, double inLon2, double *displacement);
int find_min_dx(struct latLon *ll_arr, int ll_arr_point, int ll_arr_size, int *contour_arr, int contour_arr_size, int *min_dx_index, double *min_disp);
int contained(int nps, struct index_pair pair_in_question, struct index_pair *array_of_pairs_in_question);

int wt(char i_filename[]);
int rt(char i_filename[]);
int vt(char i_filename[]);
int mrt(char i_filename[]);
int hrt(char i_filename[]);

void Contour(float **d,int ilb,int iub,int jlb,int jub,
             int *x,int *y,int nc,float *z, struct index_pair ***indp1, struct index_pair ***indp2, struct index_pair ***indp3, struct index_pair ***indp4);
void handle_error(int errid);
void malloc_error(char *In);
