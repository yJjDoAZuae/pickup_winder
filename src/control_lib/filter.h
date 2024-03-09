#ifdef __arm__
    #include "pico/stdlib.h"
#else
    #include <stdlib.h>
    #include <stdbool.h>
#endif

#define ARMA_SISO_MAX_ORDER 8
#define FIR_SISO_MAX_ORDER 128

#include "math_constants.h"

#ifndef FILTER_H
#define FILTER_H

typedef float arma_buffer_t[ARMA_SISO_MAX_ORDER+1];

// Init a buffer to zeros
int arma_buffer_init(arma_buffer_t *buff);

typedef struct
{

    int n;  // filter order (0th order is a static gain defined by b[0])
    arma_buffer_t a; // denominator coefficients, ordered by delay (z^-1).  E.g. a[0] is zero delay.  Note a[0] is ignored and set to 1.0f in all cases.
    arma_buffer_t b; // numerator coefficients, ordered by delay (z^-1). E.g. b[0] is zero delay.
    arma_buffer_t u; // input buffer
    arma_buffer_t y; // output buffer

} arma_siso_filter_state_t;

// Initialize a filter to zero state and output buffers and a zero order DC passthrough numerator and denominator
int arma_siso_filter_init(arma_siso_filter_state_t *state);

// NOTE: Only strictly stable filter designs should use these functions
// Filters with a pure integrator root will return early, but
// we don't do an eigenvalue test.

// Return DC gain of a filter design
int arma_siso_filter_dcgain(float * dcgain, arma_siso_filter_state_t *state);

// First order Tustin (bilinear) discrete time IIR filter realization
int tustin_1(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz);

// Second order Tustin (bilinear) discrete time IIR filter realization 
int tustin_2(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz);

// Update the coefficients of an ARMA SISO discrete IIR filter
int arma_siso_filter_coef_update(arma_buffer_t * num, arma_buffer_t * den, int n, arma_siso_filter_state_t * state);

int arma_siso_lowpass_1_design(float dt, float tau, arma_siso_filter_state_t *state);
int arma_siso_lowpass_2_design(float dt, float wn, float zeta, float tau_zero, arma_siso_filter_state_t *state);
int arma_siso_notch_2_design(float dt, float wn, float zeta, float tau_zero, arma_siso_filter_state_t *state);
int arma_siso_highpass_1_design(float dt, float tau, arma_siso_filter_state_t *state);
int arma_siso_highpass_2_design(float dt, float wn, float zeta, float tau_zero, arma_siso_filter_state_t *state);
int arma_siso_deriv_1_design(float dt, float tau, arma_siso_filter_state_t *state);

// Reset a filter based on a constant input value
int arma_siso_filter_input_reset(float in, arma_siso_filter_state_t *state);

// Reset a filter based on a constant output value, assuming DC gain is non-zero
int arma_siso_filter_output_reset(float out, arma_siso_filter_state_t *state);

// Step a filter forward one iteration
int arma_siso_filter_update(float in, float *out, arma_siso_filter_state_t *state);


typedef int fir_value_t;

typedef fir_value_t fir_buffer_t[FIR_SISO_MAX_ORDER];

// Init a buffer to zeros
int fir_buffer_init(fir_buffer_t *buff);

typedef struct
{

    int n;  // filter order (0th order is a static gain defined by b[0])
    fir_buffer_t b; // numerator coefficients, ordered by delay (z^-1). E.g. b[0] is zero delay.
    fir_buffer_t u; // input buffer

} fir_siso_filter_state_t;

// Return DC gain of a FIR filter design
int fir_siso_filter_dcgain(fir_value_t * dcgain, fir_siso_filter_state_t *state);

// Initialize a FIR filter
int fir_siso_filter_init(fir_siso_filter_state_t *state);

// Reset a filter based on a constant input value
int fir_siso_filter_input_reset(fir_value_t in, fir_siso_filter_state_t *state);

// Reset a filter based on an output value
int fir_siso_filter_output_reset(fir_value_t out, fir_siso_filter_state_t *state);

// Step a filter forward one iteration
int fir_siso_filter_update(fir_value_t in, fir_value_t *out, fir_siso_filter_state_t *state);

#endif
