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
Quality assurance header
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#ifndef QUALITY_CL_H
#define QUALITY_CL_H

#include <stdio.h>   // core input and output functions
#include <stdlib.h>  // standard general utilities library

#include "../cross-level/const-cl.h"
#include "../cross-level/brick_base-cl.h"


#ifdef __cplusplus
extern "C" {
#endif

/** This function sets any quality bit in the QAI layer
+++ Attention: this function implements no safety measures!
+++ The brick short memory, bitfield, pixel, and index are all vulnerable
+++ to misuse. Take care.
--- qai:    Quality Assurance Information
--- index:  QAI layer
--- p:      pixel
--- val:    set to this value (typically 0 or 1, but can be another
            integer, too, in which case, a wider bit field is changed
--- bitfields: how many bitfields to set? (typically 1 for binary bit,
               but can be larger to set multi-bit flags)
+++ Return: void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
static inline void set_qai(brick_t *qai, int index, int p, short val, int bitfields){

    // Clear the existing bits for the specified bitfields
    short mask = (1 << bitfields) - 1;
    qai->vshort[0][p] &= ~(mask << index);

    // Set the new value
    qai->vshort[0][p] |= (short)(val << index);
}

char get_aerosol(brick_t *qai, int p);
bool get_subzero(brick_t *qai, int p);
bool get_saturation(brick_t *qai, int p);
bool get_lowsun(brick_t *qai, int p);
bool get_slope(brick_t *qai, int p);
bool get_vaporfill(brick_t *qai, int p);

bool get_off_from_value(short qai_value);
char get_cloud_from_value(short qai_value);
bool get_shadow_from_value(short qai_value);
bool get_snow_from_value(short qai_value);
bool get_water_from_value(short qai_value);
char get_aerosol_from_value(short qai_value);
bool get_subzero_from_value(short qai_value);
bool get_saturation_from_value(short qai_value);
bool get_lowsun_from_value(short qai_value);
char get_illumination_from_value(short qai_value);
bool get_slope_from_value(short qai_value);
bool get_vaporfill_from_value(short qai_value);

void set_off(brick_t *qai, int p, short val);

void set_snow(brick_t *qai, int p, short val);
void set_water(brick_t *qai, int p, short val);
void set_aerosol(brick_t *qai, int p, short val);
void set_subzero(brick_t *qai, int p, short val);
void set_saturation(brick_t *qai, int p, short val);
void set_lowsun(brick_t *qai, int p, short val);
void set_illumination(brick_t *qai, int p, short val);
void set_slope(brick_t *qai, int p, short val);
void set_vaporfill(brick_t *qai, int p, short val);

void set_off_to_value(short *value, short val);
void set_shadow_to_value(short *value, short val);
void set_snow_to_value(short *value, short val);
void set_water_to_value(short *value, short val);
void set_aerosol_to_value(short *value, short val);
void set_subzero_to_value(short *value, short val);
void set_saturation_to_value(short *value, short val);
void set_lowsun_to_value(short *value, short val);
void set_illumination_to_value(short *value, short val);
void set_slope_to_value(short *value, short val);
void set_vaporfill_to_value(short *value, short val);

void merge_qai_from_values(brick_t *qai, int p, short qai_1, short qai_2);
/** This function reads any quality bit in the QAI layer
+++ Attention: this function implements no safety measures!
+++ The brick short memory, bitfield, pixel, and index are all vulnerable
+++ to misuse. Take care.
--- qai:       Quality Assurance Information
--- index:     QAI layer
--- p:         pixel
--- bitfields: how many bitfields to read? (typically 1 for binary bit,
               but can be larger to retrieve multi-bit flags
+++ Return:    void
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
static inline short get_qai(brick_t *qai, int index, int p, int bitfields){
    int i;
    short val = 0;

    for (i=0; i<bitfields; i++) val |= (short)(1 << i);

    return (short)(qai->vshort[0][p] >> index) & val;
}
/** This function sets any quality bit in the QAI layer
+++ The same as set_qai, but writes the value to a short value directly
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
static inline void set_qai_to_value(short *value, int index, short val, int bitfields){

    // Clear the existing bits for the specified bitfields
    short mask = (1 << bitfields) - 1;
    *value &= ~(mask << index);

    // Set the new value
    *value |= (short)(val << index);
}
/** This function reads any quality bit in the QAI layer
+++ The same as get_qai, but reads the value from a short value directly
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
static inline short get_qai_from_value(short value, int index, int bitfields){
    int i;
    short val = 0;

    for (i=0; i<bitfields; i++) val |= (short)(1 << i);

    return (short)(value >> index) & val;
}
/** read illumination flag **/
static inline char get_illumination(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_ILL_, p, 2);
}
/** read water flag **/
static inline bool get_water(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_WTR_, p, 1);
}
/** read snow flag **/
static inline bool get_snow(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_SNW_, p, 1);
}
/** read cloud shadow flag **/
static inline bool get_shadow(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_SHD_, p, 1);
}
/** set cloud shadow flag **/
static inline void set_shadow(brick_t *qai, int p, short val){
    set_qai(qai, _QAI_BIT_SHD_, p, val, 1);
}
/** read cloud flag **/
static inline char get_cloud(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_CLD_, p, 2);
}
/** set cloud flag, directly to value **/
static inline void set_cloud_to_value(short *value, short val){
    set_qai_to_value(value, _QAI_BIT_CLD_, val, 2);
}
/** set cloud flag **/
static inline void set_cloud(brick_t *qai, int p, short val){
    set_qai(qai, _QAI_BIT_CLD_, p, val, 2);
}
/** read off/on flag **/
static inline bool get_off(brick_t *qai, int p){

    return get_qai(qai, _QAI_BIT_OFF_, p, 1);
}

#ifdef __cplusplus
}
#endif

#endif

