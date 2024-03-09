

#ifdef __arm__
    #include "pico/stdlib.h"
#else
    #include <stdlib.h>
    #include <stdbool.h>
#endif

#define MAX_AW_STATES 12


#include "math_constants.h"
#include "filter.h"

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
