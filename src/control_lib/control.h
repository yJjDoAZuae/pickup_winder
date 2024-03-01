

#ifdef __arm__
    #include "pico/stdlib.h"
#else
    #include <stdlib.h>
    #include <stdbool.h>
#endif

#define ARMA_SISO_MAX_ORDER 8
#define MAX_AW_STATES 12
#define FIR_SISO_MAX_ORDER 128


#include "math_constants.h"

#ifndef CONTROL_H
#define CONTROL_H

// float version of clip
float fclip(float in, float vmin, float vmax);

// double version of clip
double dclip(double in, double vmin, double vmax);

// wraps an angle or angle difference,
// maps (-inf, inf) onto [-pi, pi)
float fpiwrap(float in);

// wraps an angle or angle difference,
// maps (-inf, inf) onto [-pi, pi)
double dpiwrap(double in);

// unwraps an angle based on an assumption that the previous 
// unwrapped value, when wrapped should be within a [-pi,pi) 
// interval of the current input,
// maps [-pi,pi) x (-inf,inf) onto (-inf, inf)
float fpiunwrap(float in, float prev);

// unwraps an angle based on an assumption that the previous 
// unwrapped value, when wrapped should be within a [-pi,pi) 
// interval of the current input,
// maps [-pi,pi) x (-inf,inf) onto (-inf, inf)
double dpiunwrap(double in, double prev);

typedef struct
{
    float k[ARMA_SISO_MAX_ORDER+1];
} arma_buffer_t;

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

// Step a filter forward one iteration
int arma_siso_filter_update(float in, float *out, arma_siso_filter_state_t *state);


typedef int fir_value_t;

typedef struct
{
    fir_value_t k[FIR_SISO_MAX_ORDER];
} fir_buffer_t;

// Init a buffer to zeros
int fir_buffer_init(fir_buffer_t *buff);

typedef struct
{

    int n;  // filter order (0th order is a static gain defined by b[0])
    fir_buffer_t b; // numerator coefficients, ordered by delay (z^-1). E.g. b[0] is zero delay.
    fir_buffer_t u; // input buffer

} fir_siso_filter_state_t;


typedef struct
{
    float dt;

    float min;
    float max;

    bool lower;
    bool upper;

    float y;

} limit_state_t;

int limit_init(limit_state_t *state);

int val_limit_update(float in, float * out, limit_state_t * state);
int rate_limit_reset(float in, limit_state_t *state);
int rate_limit_update(float in, float * out, limit_state_t * state);

typedef struct
{
    arma_siso_filter_state_t filter;
    limit_state_t rlim;
    limit_state_t vlim;
} signal_state_t;

int signal_init(signal_state_t *state);

int signal_reset(float in, signal_state_t *state);
int signal_update(float in, float * out, signal_state_t * state);

typedef enum {
    AW_LOWER = -1u,
    AW_ANY = 0u,
    AW_UPPER = 1u
} antiwindup_direction_t;

typedef struct
{

    // determines which direction the antiwindup mode controls
    antiwindup_direction_t dir[MAX_AW_STATES];

    // signifies whether the antiwindup mode is locked
    bool k[MAX_AW_STATES];

} antiwindup_state_t;

int antiwindup_init(antiwindup_state_t *state);

typedef struct
{

    antiwindup_state_t aw;
    limit_state_t vlim;

    float dt;

} integrator_state_t;

int integrator_init(integrator_state_t *state);

int integrator_reset(float val, integrator_state_t *state);
int integrator_update(float in, float * out, integrator_state_t * state);

typedef float gain_state_t;

typedef struct
{

    signal_state_t cmd;
    signal_state_t meas;
    signal_state_t err;
    signal_state_t out;
    signal_state_t meas_dot;

    integrator_state_t integ;

    gain_state_t Kff;
    gain_state_t Kp;
    gain_state_t Ki;
    gain_state_t Kd;

    // approximate derivative filter for generating measurement derivative
    // in cases where there is not an independent derivative measurement
    arma_siso_filter_state_t meas_deriv_filt;

    float deriv_tau;

    float dt;

} PID_state_t;

int PID_init(PID_state_t * state);

int PID_set_dt(float dt, PID_state_t * state);

int PID_reset(float cmd, float meas, float out, float integ, PID_state_t * state);
int PID_update_deriv(float cmd, float meas, float meas_dot, float * out, PID_state_t * state);
int PID_update(float cmd, float meas, float * out, PID_state_t * state);

#endif // CONTROL_H
