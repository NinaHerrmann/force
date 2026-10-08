/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This file is part of FORCE - Framework for Operational Radiometric 
Correction for Environmental monitoring.

Copyright (C) 2013-2022 David Frantz

FORCE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FORCE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FORCE.  If not, see <http://www.gnu.org/licenses/>.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions for handling topographic effects
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "topo-ll.h"


top_t *allocate_topography();
int init_topography(top_t *top);
int ocean_topography(brick_t *DEM);
int smooth_topography(brick_t *DEM);
int exposition_topography(brick_t *DEM, brick_t *EXP, brick_t *QAI);
int stats_topography(atc_t *atc, brick_t *DEM, brick_t *CDEM, brick_t *QAI);
int illumination_topography(atc_t *atc, brick_t *EXP, brick_t *ILL, brick_t *SKY, brick_t *QAI);


/** This function allocates the topographic variables
+++ Return: topographic variables (must be freed with free_topography)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
top_t *allocate_topography(){
top_t *top = NULL;


  alloc((void**)&top, 1, sizeof(top_t));
  init_topography(top);    

  return top;
}


/** This function initializes the topographic variables
--- top:    topographic variables
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int init_topography(top_t *top){

  top->dem = NULL;
  top->exp = NULL;
  top->ill = NULL;
  top->sky = NULL;
  top->c   = NULL;

  return SUCCESS;
}


/** This function detects oceans and sets these pixels to 0m a.s.l. as 
+++ ocean is often masked out in DEMs (but is valid data)
--- DEM:    Digital Elevation Model
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int ocean_topography(brick_t *DEM){
int i, j, ii, jj, p, np, nx, ny;
float *dem_ = NULL;
float nodata;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif
  
  
  nx = get_brick_ncols(DEM);
  ny = get_brick_nrows(DEM);
  nodata = get_brick_nodata(DEM, 0);
  if ((dem_ = get_band_float(DEM, 0)) == NULL) return FAILURE;

  #pragma omp parallel private(i, ii, p, np) shared(nx, ny, dem_, nodata) default(none) 
  {

    #pragma omp for schedule(static)
    for (j=0; j<nx; j++){
    for (i=0, ii=1; ii<ny; i++, ii++){
      p = i*nx+j; np = ii*nx+j;
      if (fequal(dem_[p], 0) && fequal(dem_[np], nodata)) dem_[np] = 0.0;
    }
    }

  }
  

  #pragma omp parallel private(i, ii, p, np) shared(nx, ny, dem_, nodata) default(none) 
  {

    #pragma omp for schedule(static)
    for (j=0; j<nx; j++){
    for (i=(ny-1), ii=(ny-2); ii>=0; i--, ii--){
      p = i*nx+j; np = ii*nx+j;
      if (fequal(dem_[p], 0) && fequal(dem_[np], nodata)) dem_[np] = 0.0;
    }
    }
  
  }

  
  #pragma omp parallel private(j, jj, p, np) shared(nx, ny, dem_, nodata) default(none) 
  {

    #pragma omp for schedule(static)
    for (i=0; i<ny; i++){
    for (j=0, jj=1; jj<nx; j++, jj++){
      p = i*nx+j; np = i*nx+jj;
      if (fequal(dem_[p], 0) && fequal(dem_[np], nodata)) dem_[np] = 0.0;
    }
    }
  
  }

  
  #pragma omp parallel private(j, jj, p, np) shared(nx, ny, dem_, nodata) default(none) 
  {

    #pragma omp for schedule(static)
    for (i=0; i<ny; i++){
    for (j=(nx-1), jj=(nx-2); jj>=0; j--, jj--){
      p = i*nx+j; np = i*nx+jj;
      if (fequal(dem_[p], 0) && fequal(dem_[np], nodata)) dem_[np] = 0.0;
    }
    }
    
  }


  #ifdef FORCE_CLOCK
  proctime_print("detect ocean from topography", TIME);
  #endif
  
  return SUCCESS;
}


/** This function smooths the DEM with a lowpass filter as there are often
+++ undesired effects in high resolution DEMs
--- DEM:    Digital Elevation Model
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int smooth_topography(brick_t *DEM){
int i, j, ii, jj, ni, nj, p, np, nx, ny, nc, k;
float *buf = NULL;
float sum, num;
double res;
float *dem_ = NULL;
float nodata;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif
  
  
  nx = get_brick_ncols(DEM);
  ny = get_brick_nrows(DEM);
  nc = get_brick_ncells(DEM);
  res = get_brick_res(DEM);
  nodata = get_brick_nodata(DEM, 0);
  if ((dem_ = get_band_float(DEM, 0)) == NULL) return FAILURE;


  if (res >= 30) return SUCCESS;

  #ifdef FORCE_DEV
  printf("smooth_topo should consider actual pixel size of DEM...\n");
  #endif

  k = 30.0 / res;

  alloc((void**)&buf, nc, sizeof(float));

  #pragma omp parallel private(j, p, ii, jj, ni, nj, np, sum, num) shared(nx, ny, k, dem_, buf, nodata) default(none) 
  {

    #pragma omp for schedule(guided)
    for (i=0; i<ny; i++){
    for (j=0; j<nx; j++){

      p = i*nx+j;

      buf[p] = nodata;

      if (dem_[p] == nodata) continue;

      sum = num = 0;

      for (ii=-1*k; ii<=k; ii++){
      for (jj=-1*k; jj<=k; jj++){

        ni = i+ii; nj = j+jj;
        if (ni > ny-1 || ni < 0 || nj > nx-1 || nj < 0) continue;

        np = ni*nx+nj;
        if (fequal(dem_[np], nodata)) continue;

        sum += dem_[np];
        num++;

      }
      }

      if (num > 0) buf[p] = sum/num;

    }
    }

  }
  
  memmove(dem_, buf, nc*sizeof(float));
  free((void*)buf);
  
  
  #ifdef FORCE_CLOCK
  proctime_print("smooth topography", TIME);
  #endif

  return SUCCESS;
}


/** This function computes slope and aspect with the Horn (1981) method
--- DEM:    Digital Elevation Model
--- EXP:    Exposition
--- QAI:    Quality Assurance Information
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int exposition_topography(brick_t *DEM, brick_t *EXP, brick_t *QAI){
int i, j, ii, jj, p, nx, ny, nc;
float devx, devy;
float tmp, dem[3][3];
ushort *slp_  = NULL;
ushort *asp_  = NULL;
float  *dem_  = NULL;
float nodata;
float slope, aspect;
double res;
bool valid;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif
  
  
  nx = get_brick_ncols(DEM);
  ny = get_brick_nrows(DEM);
  nc = get_brick_ncells(DEM);
  res = get_brick_res(DEM);
  nodata = (short)get_brick_nodata(DEM, 0);
  if ((dem_ = get_band_float(DEM, 0)) == NULL) return FAILURE;

  if ((slp_ = get_band_ushort(EXP, ZEN)) == NULL) return FAILURE;
  if ((asp_ = get_band_ushort(EXP, AZI)) == NULL) return FAILURE;


  #pragma omp parallel private(j, p, ii, jj, tmp, valid, devx, devy, slope, aspect, dem) shared(nx, ny, res, QAI, dem_, slp_, asp_, nodata) default(none) 
  {

    #pragma omp for schedule(guided)
    for (i=1; i<(ny-1); i++){
    for (j=1; j<(nx-1); j++){

      p = i*nx+j;

      if (get_off(QAI, p)) continue;

      valid = true;

      for (ii=-1; ii<=1; ii++){
      for (jj=-1; jj<=1; jj++){
        tmp = dem_[(i+ii)*nx+(j+jj)];
        if (fequal(tmp, nodata)){
          valid = false;
        } else {
          dem[ii+1][jj+1] = tmp;
        }
      }
      }


      if (valid){

        devx = ((dem[0][2]+2*dem[1][2]+dem[2][2]) - 
                (dem[0][0]+2*dem[1][0]+dem[2][0])) / (8*res);
        devy = ((dem[2][0]+2*dem[2][1]+dem[2][2]) - 
                (dem[0][0]+2*dem[0][1]+dem[0][2])) / (8*res);

        slope = atan(sqrt(devx*devx+devy*devy))*_R2D_CONV_;
        if (slope > 2) set_slope(QAI, p, true);

        if (slope > 0){
          aspect = atan2(devy, -1*devx)*_R2D_CONV_;
          if (aspect < 0){
            aspect = 90-aspect;
          } else if (aspect > 90){
            aspect = 360-aspect+90;
          } else {
            aspect = 90-aspect;
          }
        } else {
          slope = aspect = 0;
        }

        slp_[p] = (ushort)(slope*_D2R_CONV_*10000);
        asp_[p] = (ushort)(aspect*_D2R_CONV_*10000);

      } else {

        set_off(QAI, p, true);

      }

    }
    }
    
  }


  #pragma omp parallel shared(nx, nc, QAI, slp_, asp_) default(none) 
  {

    #pragma omp for schedule(static)
    for (p=0; p<nc; p+=nx){

      if (get_off(QAI, p)) continue;
      
      if (get_off(QAI, p+1)){
        set_off(QAI, p, true);
      } else {
        slp_[p] = slp_[p+1];
        asp_[p] = asp_[p+1];
      }
    }
    
  }
  
  
  #pragma omp parallel shared(nx, nc, QAI, slp_, asp_) default(none) 
  {

    #pragma omp for schedule(static)
    for (p=nx-1; p<nc; p+=nx){

      if (get_off(QAI, p)) continue;

      if (get_off(QAI, p-1)){
        set_off(QAI, p, true);
      } else {
        slp_[p] = slp_[p-1];
        asp_[p] = asp_[p-1];
      }
    }
    
  }
  
  
  #pragma omp parallel shared(nx, QAI, slp_, asp_) default(none) 
  {

    #pragma omp for schedule(static)
    for (p=0; p<nx; p++){

      if (get_off(QAI, p)) continue;

      if (get_off(QAI, p+nx)){
        set_off(QAI, p, true);
      } else {
        slp_[p] = slp_[p+nx];
        asp_[p] = asp_[p+nx];
      }
    }
    
  }
  
  
  #pragma omp parallel shared(nx, ny, nc, QAI, slp_, asp_) default(none) 
  {

    #pragma omp for schedule(static)
    for (p=(ny-1)*nx; p<nc; p++){

      if (get_off(QAI, p)) continue;

      if (get_off(QAI, p-nx)){
        set_off(QAI, p, true);
      } else {
        slp_[p] = slp_[p-nx];
        asp_[p] = asp_[p-nx];
      }
    }
    
  }

  
  #ifdef FORCE_CLOCK
  proctime_print("exposition from topography", TIME);
  #endif

  return SUCCESS;
}


/** This function computes elevation statistics (in km) and computes the
+++ binned DEM (100m classes)
--- atc:    atmospheric correction factors
--- DEM:    Digital Elevation Model
--- CDEM:   Binned Digital Elevation Model
--- QAI:    Quality Assurance Information
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int stats_topography(atc_t *atc, brick_t *DEM, brick_t *CDEM, brick_t *QAI){
int p, nc;
float dem, mn = SHRT_MAX, mx = SHRT_MIN;
double sum = 0, num = 0;
float  *dem_  = NULL;
small  *cdem_ = NULL;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif
  

  nc = get_brick_ncells(DEM);
  if ((dem_ = get_band_float(DEM, 0)) == NULL) return FAILURE;
  if ((cdem_ = get_band_small(CDEM, 0)) == NULL) return FAILURE;


  #pragma omp parallel private(dem) shared(nc, QAI, dem_) reduction(+: sum, num) reduction(max: mx) reduction(min: mn) default(none) 
  {

    #pragma omp for schedule(guided)
    for (p=0; p<nc; p++){
      if (get_off(QAI, p)) continue;
      dem = dem_[p]/1000.0; // -> kilometer
      if (dem < mn) mn = dem;
      if (dem > mx) mx = dem;
      sum += dem; // -> kilometer
      num++;
    }
    
  }

  atc->dem.min = mn;
  atc->dem.max = mx;
  atc->dem.avg = (float)(sum/num);
  atc->dem.max += 0.001; // add 1m

  atc->dem.cnum = _BYTE_LEN_ - 1;
  atc->dem.step = (atc->dem.max-atc->dem.min)/atc->dem.cnum;
  
  
  #pragma omp parallel private(dem) shared(nc, QAI, dem_, cdem_, atc) default(none) 
  {

    #pragma omp for schedule(guided)
    for (p=0; p<nc; p++){
      if (get_off(QAI, p)) continue;
      dem = dem_[p]/1000.0; // -> kilometer
      cdem_[p] = (small)floor((dem-atc->dem.min)/atc->dem.step);
    }

  }


  if (atc->dem.min < -0.5 || atc->dem.max > 9){
    printf("DEM out of bounds: min %f max %f. ", 
    atc->dem.min, atc->dem.max); return FAILURE;
  }

 
  #ifdef FORCE_DEBUG
  printf("elevation stats:\n");
  printf(" avg: %+6.3f, min: %+6.3f, max: %+6.3f, step: %+6.3f, cnum: %d\n",
  atc->dem.avg, atc->dem.min, atc->dem.max, atc->dem.step, atc->dem.cnum);
  #endif

  #ifdef FORCE_CLOCK
  proctime_print("stats topography", TIME);
  #endif

  return SUCCESS;
}


/** This function computes the illumination angle and a simple sky view 
+++ factor
--- atc:    atmospheric correction factors
--- EXP:    Exposition
--- ILL:    Illumination angle
--- SKY:    Sky view factor
--- QAI:    Quality Assurance Information
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int illumination_topography(atc_t *atc, brick_t *EXP, brick_t *ILL, brick_t *SKY, brick_t *QAI){
int p, g, nc;
float slp, asp;
ushort *slp_ = NULL;
ushort *asp_ = NULL;
short  *ill_ = NULL;
ushort *sky_ = NULL;
float **sun_ = NULL;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif


  nc = get_brick_ncells(EXP);
  if ((slp_ = get_band_ushort(EXP, ZEN)) == NULL) return FAILURE;
  if ((asp_ = get_band_ushort(EXP, AZI)) == NULL) return FAILURE;
  if ((ill_ = get_band_short(ILL, 0)) == NULL) return FAILURE;
  if ((sky_ = get_band_ushort(SKY, 0)) == NULL) return FAILURE;
  if ((sun_ = get_bands_float(atc->xy_sun)) == NULL) return FAILURE;


  #pragma omp parallel private(g, slp, asp) shared(nc, QAI, slp_, asp_, ill_, sky_, sun_, atc) default(none) 
  {

    #pragma omp for schedule(guided)
    for (p=0; p<nc; p++){

      if (get_off(QAI, p)){ ill_[p] = -10000; continue;}

      g = convert_brick_p2p(QAI, atc->xy_sun, p);

      slp = slp_[p]/10000.0;
      asp = asp_[p]/10000.0;

      // illumination angle
      if (slp == 0){
        ill_[p] = (short)(sun_[cZEN][g]*10000);
      } else {
        ill_[p] = (short)(illumin(sun_[cZEN][g], sun_[sZEN][g], 
                              cos(slp), sin(slp), sun_[AZI][g], asp)*10000);
      }

      // set illumination QAI
      if (ill_[p] < 0){
        set_illumination(QAI, p, 3); // deep shadow
      } else if (ill_[p] < 1736.482){
        set_illumination(QAI, p, 2); // poor
      } else if (ill_[p] < 5735.764){
        set_illumination(QAI, p, 1); // moderate
      }

      // sky view factor (portion of the sky dome diffusing on to a tilted surface)
      sky_[p] = (ushort)(10000 - (slp/M_PI)*10000);

    }
    
  }


  #ifdef FORCE_CLOCK
  proctime_print("illumination condition", TIME);
  #endif

  return SUCCESS;
}


/** public functions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


/** This function frees the topographic variables
--- top:    topographic variables
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void free_topography(top_t *top){

  if (top == NULL) return;

  free_brick(top->dem);
  free_brick(top->exp);
  free_brick(top->ill);
  free_brick(top->sky);
  free_brick(top->c);

  free((void*)top); top = NULL;

  return;
}


/** This function compiles the basic set of topographic derivatives. This 
+++ includes reprojecting the DEM to the image extent/projection and com-
+++ puting slope and aspect
--- pl2:        L2 parameters
--- atc:        atmospheric correction factors
--- topography: topographic variables
--- QAI:        Quality Assurance Information
+++ Return:     SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int compile_topography(par_ll_t *pl2, atc_t *atc, top_t **topography, brick_t *QAI){
top_t *top = NULL;
brick_t *DEM = NULL;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif
  
  top = allocate_topography();


  /** Digital Elevation Model brick **/
  DEM = copy_brick(QAI, 1, _DT_FLOAT_);
  set_brick_name(DEM, "FORCE DEM brick");
  set_brick_product(DEM, "DEM");
  set_brick_filename(DEM, "DEM");
  set_brick_bandname(DEM, 0, "DEM");
  set_brick_nodata(DEM, 0, pl2->dem_nodata);


  // warp DEM to MEM or use flat DEM (z = 0m)
  if (strcmp(pl2->fdem, "NULL") != 0){
    if ((warp_from_disc_to_known_brick(pl2->dem_resample, pl2->nthread, pl2->fdem, DEM, 0, 0, pl2->dem_nodata)) != SUCCESS){
      printf("Reprojecting of DEM failed! "); return FAILURE;}
  }


  /** detect oceans and set to 0m a.s.l **/
  if ((ocean_topography(DEM)) != SUCCESS){
    printf("Compiling ocean DEM failed! "); return FAILURE;}


  /** smooth DEM **/
  if ((smooth_topography(DEM)) != SUCCESS){
    printf("Smoothing of DEM failed! "); return FAILURE;}


  #ifdef FORCE_DEBUG
  print_brick_info(DEM); set_brick_open(DEM, OPEN_CREATE); write_brick(DEM);
  #endif


  /** exposition brick **/
  top->exp = copy_brick(DEM, 2, _DT_USHORT_);
  set_brick_name(top->exp, "FORCE terrain exposition brick");
  set_brick_product(top->exp, "EXP");
  set_brick_filename(top->exp, "DEM-EXPOSITION");
  set_brick_bandname(top->exp, ZEN, "Slope");
  set_brick_bandname(top->exp, AZI, "Aspect");

  if ((exposition_topography(DEM, top->exp, QAI)) != SUCCESS){
    printf("Slope/aspect failed! "); return FAILURE;}
    
  #ifdef FORCE_DEBUG
  print_brick_info(top->exp); set_brick_open(top->exp, OPEN_CREATE); write_brick(top->exp);
  #endif

  
  /** calculate DEM stats + binned DEM **/
  top->dem = copy_brick(DEM, 1, _DT_SMALL_);
  set_brick_name(top->dem, "FORCE binned DEM brick");
  set_brick_product(top->dem, "BEM");
  set_brick_filename(top->dem, "DEM-BINNED");
  set_brick_bandname(top->dem, 0, "binned DEM");

  if ((stats_topography(atc, DEM, top->dem, QAI)) != SUCCESS){
    printf("Elevation statistics failed! "); return FAILURE;}

  #ifdef FORCE_DEBUG
  print_brick_info(top->dem); set_brick_open(top->dem, OPEN_CREATE); write_brick(top->dem);
  #endif


  /** illumination angle and sky view factor **/
  top->ill = copy_brick(DEM, 1, _DT_SHORT_);
  set_brick_name(top->ill, "FORCE Illumination angle brick");
  set_brick_product(top->ill, "ILL");
  set_brick_filename(top->ill, "DEM-ILLUMINATION");
  set_brick_bandname(top->ill, 0, "Illumination angle");

  top->sky = copy_brick(DEM, 1, _DT_USHORT_);
  set_brick_name(top->sky, "FORCE Sky View Factor brick");
  set_brick_product(top->sky, "SKY");
  set_brick_filename(top->sky, "DEM-SKY-VIEW");
  set_brick_bandname(top->sky, 0, "Sky View Factor");

  if (illumination_topography(atc, top->exp, top->ill, top->sky, QAI) != SUCCESS){
    printf("error in topographic correction. "); return FAILURE;}

  free_brick(DEM);


  #ifdef FORCE_DEBUG
  print_brick_info(top->ill); set_brick_open(top->ill, OPEN_CREATE); write_brick(top->ill);
  print_brick_info(top->sky); set_brick_open(top->sky, OPEN_CREATE); write_brick(top->sky);
  print_brick_info(QAI); set_brick_open(QAI, OPEN_CREATE); write_brick(QAI);
  #endif

  #ifdef FORCE_CLOCK
  proctime_print("compiled topography", TIME);
  #endif

  *topography = top;
  return SUCCESS;
}


/** This function computes the C-factor used for topographic correction. C
+++ is derived for the SWIR2 band and is then propagated through the spec-
+++ trum in the actual radiometric correction. C is derived for each pixel
+++ with sufficient illumination.
--- atc:    atmospheric correction factors
--- TOA:    TOA reflectance
--- QAI:    Quality Assurance Information
--- DEM:    Digital Elevation model
--- EXP:    Exposition
--- ILL:    Illumination angle
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
brick_t *cfactor_topography(atc_t *atc, brick_t *TOA, brick_t *QAI, brick_t *DEM, brick_t *EXP, brick_t *ILL){
int i, j, p, ip, jp, np, nx, ny, nc;
int b_sw2;
ushort s_min  = 350; // 2° slope
double mx;
double  num;
float cf;
float *cor_ = NULL;
float tmp;
double res;
float *swir_ = NULL;
int k, nk = 0, *K = NULL;
brick_t *CF = NULL;
ushort  *cf_ = NULL;
short *sw1_ = NULL;
short *sw2_ = NULL;
small *dem_ = NULL;
ushort *slp_ = NULL;
short *ill_ = NULL;
float *xy_szen = NULL;
float *xy_ms = NULL;
float **xyz_rho_p = NULL;
float **xyz_tss = NULL;
float **xyz_tsd = NULL;
double *xd = NULL, *yd = NULL, *wtab = NULL;
float  *swn = NULL;


  #ifdef FORCE_CLOCK
  time_t TIME; time(&TIME);
  #endif


  cite_me(_CITE_TOPCOR_);


  CF = copy_brick(QAI, 1, _DT_USHORT_);
  set_brick_name(CF, "FORCE C-factor brick");
  set_brick_product(CF, "CFC");
  set_brick_filename(CF, "DEM-C-FACTOR");
  set_brick_bandname(CF, 0, "C-Factor");
  set_brick_nodata(CF, 0, _FORCE_NO_DATA_);

  nx  = get_brick_ncols(QAI);
  ny  = get_brick_nrows(QAI);
  nc  = get_brick_ncells(QAI);
  res = get_brick_res(QAI);

  if ((cf_  = get_band_ushort(CF, 0))         == NULL) return NULL;
  if ((sw1_ = get_domain_short(TOA, "SWIR1")) == NULL) return NULL;
  if ((sw2_ = get_domain_short(TOA, "SWIR2")) == NULL) return NULL;
  if ((dem_ = get_band_small(DEM, 0))         == NULL) return NULL;
  if ((slp_ = get_band_ushort(EXP, ZEN))      == NULL) return NULL;
  if ((ill_ = get_band_short(ILL, 0))         == NULL) return NULL;
  if ((b_sw2 = find_domain(TOA, "SWIR2")) < 0) return NULL;

  if ((xy_szen     = get_band_float(atc->xy_sun,   ZEN)) == NULL) return NULL;
  if ((xy_ms       = get_band_float(atc->xy_sun,  cZEN)) == NULL) return NULL;
  if ((xyz_rho_p       = atc_get_band_reshaped(atc->xyz_rho_p, b_sw2))     == NULL) return NULL;
  if ((xyz_tss       = atc_get_band_reshaped(atc->xyz_tss, b_sw2))     == NULL) return NULL;
  if ((xyz_tsd       = atc_get_band_reshaped(atc->xyz_tsd, b_sw2))     == NULL) return NULL;


  /** kernel for sampling neighborhood **/
  alloc((void**)&K, 3000/res*2+1, sizeof(int));

  K[nk] = 3000/res; nk++;
  while (K[nk-1] > 100/res){
    K[nk] = K[nk-1]/sqrt(2);
    nk++;
  }

  if (K[nk-1] != 0) nk++;

  for (k=nk-2; k>=0; k--){
    K[nk] = -1*K[k];
    nk++;
  }

  #ifdef FORCE_DEBUG
  printf("sampling for estimating C:\n");
  for (k=0; k<nk; k++) printf("%d ", K[k]); 
  printf("\n");
  #endif


  /** allocate memory **/
  alloc((void**)&cor_,  nc, sizeof(float));
  alloc((void**)&swir_, nc, sizeof(float));
  alloc((void**)&xd,   nc, sizeof(double));
  alloc((void**)&yd,   nc, sizeof(double));
  alloc((void**)&swn,  nc, sizeof(float));
  alloc((void**)&wtab, nk*nk+1, sizeof(double));


  /** compute SWIR index **/

  #pragma omp parallel private(tmp) shared(nc, QAI, sw1_, sw2_, swir_) default(none)
  {

    #pragma omp for schedule(guided)

    for (p=0; p<nc; p++){

      if (get_off(QAI, p)) continue;

      tmp = sw1_[p]/10000.0+sw2_[p]/10000.0;
      if (tmp == 0) swir_[p] = 0; else swir_[p] = (sw1_[p]/10000.0-sw2_[p]/10000.0)/tmp;

    }

  }

  /** precompute x/y as double, wtab, and NaN-marked SWIR **/
  for (k=1; k<=nk*nk; k++) wtab[k] = (k-1)/(double)k;

  #pragma omp parallel for schedule(guided)
  for (p=0; p<nc; p++){
      xd[p] = ill_[p]/10000.0;
      yd[p] = sw2_[p]/10000.0;
      swn[p] = (!get_off(QAI, p) && ill_[p] >= 0 && slp_[p] >= s_min) ? swir_[p] : NAN;
  }
  /** estimate C for every pixel **/

  #pragma omp parallel
{
  float  *swp = NULL, *c_row = NULL, *rho_row = NULL;
  double *st  = NULL;                       // st[0..5][nx]: mx,my,vx,vy,cv,num
  alloc((void**)&swp,     nx+4, sizeof(float));
  alloc((void**)&c_row,   nx,   sizeof(float));
  alloc((void**)&rho_row, nx,   sizeof(float));
  alloc((void**)&st,      6*nx, sizeof(double));
  const int K0 = K[0];

  #pragma omp for schedule(dynamic,1)
  for (int i=0; i<ny; i++){

    /* A: per-pixel setup (same formulas as before) */
    for (int j=0; j<nx; j++){
      int p = i*nx+j;
      swp[j] = NAN; st[5*nx+j] = 0;
      if (get_off(QAI, p) || ill_[p] < 0) continue;
      int g = convert_brick_ji2p(QAI, atc->xy_sun, i, j);
      int z = dem_[p];
      float szen = xy_szen[g], ms = xy_ms[g];
      float f  = f_factor(xyz_tss[z][g], xyz_tsd[z][g]);
      float h0 = (M_PI+2*szen)/(2.0*M_PI);
      c_row[j]   = c_factor_com(h0, f, ms);
      rho_row[j] = xyz_rho_p[z][g];
      if (slp_[p] > s_min) swp[j] = swir_[p];     // NaN = no sampling
    }

    /* B: sampling */
    const int row_in = (i >= K0 && i <= ny-1-K0);
    int j = 0;
    while (j < nx){
      int p = i*nx+j;

      if (row_in && j >= K0 && j+3 <= nx-1-K0){   /* 4 pixels at once */
        const __m128 sgn = _mm_set1_ps(-0.0f);
        const __m256d thr = _mm256_set1_pd(0.025), one = _mm256_set1_pd(1.0);
        const __m128 sp = _mm_loadu_ps(swp+j);
        __m256d mx = _mm256_setzero_pd(), my = mx, vx = mx, vy = mx, cv = mx, cnt = mx;

        for (int ii=0; ii<nk; ii++){
          const long r = p + (long)K[ii]*nx;
          for (int jj=0; jj<nk; jj++){
            const long q = r + K[jj];
            __m128 d = _mm_andnot_ps(sgn, _mm_sub_ps(sp, _mm_loadu_ps(swn+q)));
            __m256d m = _mm256_cmp_pd(_mm256_cvtps_pd(d), thr, _CMP_LE_OQ);
            if (_mm256_testz_pd(m, m)) continue;  // NaN lanes never pass

            __m256d x = _mm256_loadu_pd(xd+q), y = _mm256_loadu_pd(yd+q);
            __m256d n  = _mm256_add_pd(cnt, one);
            __m256d dx = _mm256_sub_pd(x, mx), dy = _mm256_sub_pd(y, my);
            __m256d nmx = _mm256_add_pd(mx, _mm256_div_pd(dx, n));
            __m256d nmy = _mm256_add_pd(my, _mm256_div_pd(dy, n));
            __m256d nvx = _mm256_add_pd(vx, _mm256_mul_pd(dx, _mm256_sub_pd(x, nmx)));
            __m256d nvy = _mm256_add_pd(vy, _mm256_mul_pd(dy, _mm256_sub_pd(y, nmy)));
            __m256d w   = _mm256_div_pd(cnt, n);              // (n-1)/n
            __m256d ncv = _mm256_add_pd(cv, _mm256_mul_pd(_mm256_mul_pd(w, dx), dy));

            mx = _mm256_blendv_pd(mx, nmx, m);  my = _mm256_blendv_pd(my, nmy, m);
            vx = _mm256_blendv_pd(vx, nvx, m);  vy = _mm256_blendv_pd(vy, nvy, m);
            cv = _mm256_blendv_pd(cv, ncv, m);  cnt = _mm256_blendv_pd(cnt, n, m);
          }
        }
        _mm256_storeu_pd(st+0*nx+j, mx); _mm256_storeu_pd(st+1*nx+j, my);
        _mm256_storeu_pd(st+2*nx+j, vx); _mm256_storeu_pd(st+3*nx+j, vy);
        _mm256_storeu_pd(st+4*nx+j, cv); _mm256_storeu_pd(st+5*nx+j, cnt);
        j += 4;

      } else {                                     /* scalar (edges) */
        if (!isnan(swp[j])){
          double smx=0, smy=0, svx=0, svy=0, scv=0; int n = 0;
          const float sw = swp[j];
          for (int ii=0; ii<nk; ii++){
            if (i+K[ii] < 0 || i+K[ii] > ny-1) continue;
            for (int jj=0; jj<nk; jj++){
              if (j+K[jj] < 0 || j+K[jj] > nx-1) continue;
              long q = p + (long)K[ii]*nx + K[jj];
              if (!(fabs(sw-swn[q]) <= 0.025)) continue;   // NaN rejects
              double x = xd[q], y = yd[q];
              n++;
              if (n == 1){ smx = x; smy = y; }
              else {
                double dx = x-smx, dy = y-smy;
                double nmx = smx + dx/n, nmy = smy + dy/n;
                svx = svx + dx*(x-nmx);
                svy = svy + dy*(y-nmy);
                scv = scv + wtab[n]*dx*dy;
                smx = nmx; smy = nmy;
              }
            }
          }
          st[0*nx+j]=smx; st[1*nx+j]=smy; st[2*nx+j]=svx;
          st[3*nx+j]=svy; st[4*nx+j]=scv; st[5*nx+j]=n;
        }
        j++;
      }
    }

    /* C: finalize (original logic) */
    for (int j=0; j<nx; j++){
      int p = i*nx+j;
      if (get_off(QAI, p) || ill_[p] < 0) continue;
      double num = st[5*nx+j], gain, offset;
      float c_ = c_row[j];
      if (num > 2){
        double mx = st[0*nx+j], my = st[1*nx+j];
        linreg_coefs(mx, my, covariance(st[4*nx+j], num),
                     variance(st[2*nx+j], num), &gain, &offset);
        if (offset < rho_row[j]){ offset = rho_row[j]; gain = (my-offset)/mx; }
        if (offset < 10*gain) cor_[p] = c_factor_emp(offset, gain);
        else                  cor_[p] = 10.0;
      }
      if (cor_[p] < 0 || cor_[p] < c_) cor_[p] = c_;
      if (cor_[p] > USHRT_MAX/10000.0) cor_[p] = USHRT_MAX/10000.0;
    }
  }
  free((void*)swp); free((void*)c_row); free((void*)rho_row); free((void*)st);
}

  free((void*)swir_);
  free((void*)xd); free((void*)yd); free((void*)swn); free((void*)wtab);
  free((void*)xyz_rho_p); free((void*)xyz_tss); free((void*)xyz_tsd);


  /** smooth C-factor with lowpass **/

  #pragma omp parallel private(j, p, ip, jp, np, num, mx, cf) shared(nx, ny, QAI, ill_, cor_, cf_) default(none)
  {

    #pragma omp for schedule(guided)
    for (i=0; i<ny; i++){
    for (j=0; j<nx; j++){

      p = i*nx+j;

      if (get_off(QAI, p) || ill_[p] < 0) continue;

        num = mx = 0.0;

        for (ip=(i-1); ip<=(i+1); ip++){
        for (jp=(j-1); jp<=(j+1); jp++){

          if (ip < 0 || jp < 0 || ip > ny-1 || jp > nx-1) continue;
          np = ip*nx+jp;

          if (get_off(QAI, np) || ill_[np] < 0) continue;

          mx += cor_[np];
          num++;

        }
        }

        if (num > 0) cf = mx/num*10000; else cf = cor_[p]*10000;
        if (cf > USHRT_MAX) cf = USHRT_MAX;
        if (num > 0) cf_[p] = (ushort)cf;

    }
    }

  }

  free((void*)cor_);
  free((void*)K);


  #ifdef FORCE_DEBUG
  print_brick_info(CF); set_brick_open(CF, OPEN_CREATE); write_brick(CF);
  #endif

  #ifdef FORCE_CLOCK
  proctime_print("topographic correction factors", TIME);
  #endif

  return CF;
}

/** Average elevation of coarse grid cell
+++ This function computes the average binned elevation in the given
+++ coarse grid cell
--- g:   cell of coarse grid
--- CDEM:   Binned Digital Elevation Model
--- FDEM:   Float  Digital Elevation Model
--- QAI:    Quality Assurance Information
+++ Return: SUCCESS/FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
int average_elevation_cell(int g, brick_t *CDEM, brick_t *FDEM, brick_t *QAI){
int i, j, ii, jj, p, cell_size, nx, ny;
double sum = 0, num = 0;
small *dem_ = NULL;


  // coarse cell to fine pos
  convert_brick_p2ji(CDEM, FDEM, g, &i, &j);

  // cellsize in fine pixels
  cell_size = floor(get_brick_res(CDEM)/get_brick_res(FDEM));

  nx = get_brick_ncols(FDEM);
  ny = get_brick_nrows(FDEM);

  if ((dem_  = get_band_small(FDEM, 0)) == NULL) return FAILURE;

  for (ii=i; ii<(i+cell_size); ii++){
  for (jj=j; jj<(j+cell_size); jj++){

    if (ii > ny-1 || jj > nx-1) continue;

    p = ii*nx+jj;
    if (get_off(QAI, p)) continue;

    sum += dem_[p];
    num++;

  }
  }

  if (num>0){
    set_brick(CDEM, 0, g, round(sum/num));
  } else {
    set_brick(CDEM, 0, g, 0);
  }
  
  return SUCCESS;
}

