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
This file contains functions for quality assurance
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "quality-cl.h"








/** read off/on flag, directly from value **/
bool get_off_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_OFF_, 1);
}



/** read cloud flag, directly from value **/
char get_cloud_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_CLD_, 2);
}



/** read cloud shadow flag, directly from value **/
bool get_shadow_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SHD_, 1);
}



/** read snow flag, directly from value **/
bool get_snow_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SNW_, 1);
}



/** read water flag, directly from value **/
bool get_water_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_WTR_, 1);
}

/** read aerosol flag **/
char get_aerosol(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_AOD_, p, 2);
}

/** read aerosol flag, directly from value **/
char get_aerosol_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_AOD_, 2);
}

/** read subzero reflectance flag **/
bool get_subzero(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_SUB_, p, 1);
}

/** read subzero reflectance flag, directly from value **/
bool get_subzero_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SUB_, 1);
}

/** read saturated reflectance flag **/
bool get_saturation(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_SAT_, p, 1);
}

/** read saturated reflectance flag, directly from value **/
bool get_saturation_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SAT_, 1);
}

/** read low sun angle flag **/
bool get_lowsun(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_SUN_, p, 1);
}

/** read low sun angle flag, directly from value **/
bool get_lowsun_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SUN_, 1);
}



/** read illumination flag, directly from value **/
char get_illumination_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_ILL_, 2);
}

/** read slope flag **/
bool get_slope(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_SLP_, p, 1);
}

/** read slope flag, directly from value **/
bool get_slope_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_SLP_, 1);
}

/** read water vapor fill flag **/
bool get_vaporfill(brick_t *qai, int p){

  return get_qai(qai, _QAI_BIT_WVP_, p, 1);
}

/** read water vapor fill flag, directly from value **/
bool get_vaporfill_from_value(short value){

  return get_qai_from_value(value, _QAI_BIT_WVP_, 1);
}

/** set off/on flag **/
void set_off(brick_t *qai, int p, short val){

  set_qai(qai, _QAI_BIT_OFF_, p, val, 1);
}

/** set off/on flag, directly to value **/
void set_off_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_OFF_, val, 1);
}

/** set cloud shadow flag, directly to value **/
void set_shadow_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SHD_, val, 1);
}

/** set snow flag **/
void set_snow(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_SNW_, p, val, 1);
}

/** set snow flag, directly to value **/
void set_snow_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SNW_, val, 1);
}

/** set water flag **/
void set_water(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_WTR_, p, val, 1);
}

/** set water flag, directly to value **/
void set_water_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_WTR_, val, 1);
}

/** set aerosol flag **/
void set_aerosol(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_AOD_, p, val, 2);
}

/** set aerosol flag, directly to value **/
void set_aerosol_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_AOD_, val, 2);
}

/** set subzero reflectance flag **/
void set_subzero(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_SUB_, p, val, 1);
}

/** set subzero reflectance flag, directly to value **/
void set_subzero_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SUB_, val, 1);
}

/** set saturated reflectance flag **/
void set_saturation(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_SAT_, p, val, 1);
}

/** set saturated reflectance flag, directly to value **/
void set_saturation_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SAT_, val, 1);
}

/** set low sun angle flag **/
void set_lowsun(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_SUN_, p, val, 1);
}

/** set low sun angle flag, directly to value **/
void set_lowsun_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SUN_, val, 1);
}

/** set illumination flag **/
void set_illumination(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_ILL_, p, val, 2);
}

/** set illumination flag, directly to value **/
void set_illumination_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_ILL_, val, 2);
}

/** set slope flag **/
void set_slope(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_SLP_, p, val, 1);
}

/** set slope flag, directly to value **/
void set_slope_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_SLP_, val, 1);
}

/** set water vapor fill flag **/
void set_vaporfill(brick_t *qai, int p, short val){
  set_qai(qai, _QAI_BIT_WVP_, p, val, 1);
}

/** set water vapor fill flag, directly to value **/
void set_vaporfill_to_value(short *value, short val){
  set_qai_to_value(value, _QAI_BIT_WVP_, val, 1);
}


/** This function merges two QAI values, using the most restrictive one
+++ Attention: this function implements no safety measures! 
+++ The brick short memory, bitfield, pixel, and index are all vulnerable
+++ to misuse. Take care.
--- qai:    Quality Assurance Information
--- p:      pixel
--- qai_1:  1st QAI value
--- qai_2:  2nd QAI value
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
void merge_qai_from_values(brick_t *qai, int p, short qai_1, short qai_2){

  // Merge off/on flag
  set_off(qai, p, get_off_from_value(qai_1) && get_off_from_value(qai_2));

  // Merge cloud flag
  set_cloud(qai, p, 
    (get_cloud_from_value(qai_1) > get_cloud_from_value(qai_2)) ? 
     get_cloud_from_value(qai_1) : get_cloud_from_value(qai_2));

  // Merge cloud shadow flag
  set_shadow(qai, p, 
    get_shadow_from_value(qai_1) || get_shadow_from_value(qai_2));

  // Merge snow flag
  set_snow(qai, p, 
    get_snow_from_value(qai_1) || get_snow_from_value(qai_2));

  // Merge water flag
  set_water(qai, p, 
    get_water_from_value(qai_1) || get_water_from_value(qai_2));

  // Merge aerosol flag
  set_aerosol(qai, p, 
    (get_aerosol_from_value(qai_1) > get_aerosol_from_value(qai_2)) ? 
     get_aerosol_from_value(qai_1) : get_aerosol_from_value(qai_2));

  // Merge subzero reflectance flag
  set_subzero(qai, p, 
    get_subzero_from_value(qai_1) || get_subzero_from_value(qai_2));

  // Merge saturated reflectance flag
  set_saturation(qai, p, 
    get_saturation_from_value(qai_1) || get_saturation_from_value(qai_2));

  // Merge low sun angle flag
  set_lowsun(qai, p, 
    get_lowsun_from_value(qai_1) || get_lowsun_from_value(qai_2));

  // Merge illumination flag
  set_illumination(qai, p, 
    (get_illumination_from_value(qai_1) < get_illumination_from_value(qai_2)) ? 
     get_illumination_from_value(qai_1) : get_illumination_from_value(qai_2));

  // Merge slope flag
  set_slope(qai, p, 
    get_slope_from_value(qai_1) || get_slope_from_value(qai_2));

  // Merge water vapor fill flag
  set_vaporfill(qai, p, 
    get_vaporfill_from_value(qai_1) || get_vaporfill_from_value(qai_2));

}
