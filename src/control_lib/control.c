#include "math.h"
#include <stdio.h>
#include <string.h>

#include "control.h"

float fclip(float in, float vmin, float vmax)
{
    return fminf(fmaxf(in, vmin), vmax);
}

double dclip(double in, double vmin, double vmax)
{
    return fmin(fmax(in, vmin), vmax);
}

float fpiwrap(float in)
{
    float wraps = floorf(in/TWOPI + 0.5);
    return in - TWOPI*wraps;
}

double dpiwrap(double in)
{
    double wraps = floor(in/TWOPI + 0.5);
    return in - TWOPI*wraps;
}

float fpiunwrap(float in, float prev)
{
    float del = fpiwrap(in-prev);
    return prev+del;
}

double dpiunwrap(double in, double prev)
{
    double del = dpiwrap(in-prev);
    return prev+del;
}

int signal_init(signal_state_t *state)
{

    float vlim=0.0f;

    int rc = 0;
    rc = limit_init(&(state->vlim));
    if (rc > 0) {
        return 1;
    }

    rc = limit_init(&(state->rlim));
    if (rc > 0) {
        return 2;
    }

    rc = arma_siso_filter_init(&(state->filter));
    if (rc > 0) {
        return 3;
    }

    return 0;
}

int signal_reset(float in, signal_state_t *state)
{

    float vlim=0.0f;

    int rc = 0;
    rc = val_limit_update(in, &vlim, &(state->vlim));
    if (rc > 0) {
        return 1;
    }

    rc = rate_limit_reset(vlim, &(state->rlim));
    if (rc > 0) {
        return 2;
    }

    rc = arma_siso_filter_input_reset(vlim, &(state->filter));
    if (rc > 0) {
        return 3;
    }

    return 0;
}

int signal_update(float in, float *out, signal_state_t *state)
{
    float vlim=0.0f;
    float rlim=0.0f;

    int rc = 0;
    rc = val_limit_update(in, &vlim, &(state->vlim));
    if (rc > 0) {
        return 1;
    }

    rc = rate_limit_update(vlim, &rlim, &(state->rlim));
    if (rc > 0) {
        return 2;
    }

    // NOTE: if filter is underdamped it may violate either 
    // the value limit or the rate limit
    arma_siso_filter_update(rlim, out, &(state->filter));
    if (rc > 0) {
        return 3;
    }

    return 0;
}

int antiwindup_init(antiwindup_state_t *state)
{
    for (int k = 0; k<MAX_AW_STATES; k++) {
        state->dir[k] = AW_ANY;
        state->k[k] = false;
    }

    return 0;
}

int integrator_init(integrator_state_t *state)
{
    float tmp;
    limit_init(&(state->vlim));
    antiwindup_init(&(state->aw));
    state->dt = 1.0f;

    return 0;
}

int integrator_reset(float val, integrator_state_t *state)
{
    float tmp;
    val_limit_update(val, &tmp, &(state->vlim));

    return 0;
}

int integrator_update(float in, float *out, integrator_state_t *state)
{

    if (state->dt <= 0.0f) {
        return 1;
    }

    float del = in*state->dt;

    // evaluate the antiwindup modes
    // e.g. if del is positive, state->aw.k[k]) is true, and state->aw.dir[k] is negative, 
    // then lock the integrator
    // state->aw.dir[k] defines the direction where the integrator is *unlocked*
    bool aw_locked = false;
    for (int k = 0; k<MAX_AW_STATES; k++) {
        aw_locked |= (state->aw.k[k]) && ((int)state->aw.dir[k]*in < 0);
    }

    float val = state->vlim.y;

    if (!aw_locked) {
        val += del;
    }

    val_limit_update(val, out, &(state->vlim));

    return 0;
}

int PID_init(PID_state_t *state)
{

    signal_init(&(state->cmd));
    signal_init(&(state->meas));
    signal_init(&(state->meas_dot));
    signal_init(&(state->err));
    signal_init(&(state->out));

    integrator_init(&(state->integ));

    arma_siso_filter_init(&(state->meas_deriv_filt));

    state->Kff = 0.0f;
    state->Kp = 0.0f;
    state->Ki = 0.0f;
    state->Kd = 0.0f;

    state->deriv_tau = 0.0f;
    state->dt = 1.0f;

    return 0;
}

int PID_set_dt(float dt, PID_state_t * state) {

    state->cmd.rlim.dt = dt;
    state->meas.rlim.dt = dt;
    state->meas_dot.rlim.dt = dt;
    state->err.rlim.dt = dt;
    state->out.rlim.dt = dt;
    state->integ.dt = dt;
    state->dt = dt;

    arma_siso_deriv_1_design(dt, state->deriv_tau, &(state->meas_deriv_filt));

}

int PID_reset(float cmd, float meas, float out, float integ, PID_state_t *state)
{

    signal_reset(cmd, &(state->cmd));
    signal_reset(meas, &(state->meas));
    signal_reset(0.0f, &(state->meas_dot));
    signal_reset(cmd-meas, &(state->err));
    signal_reset(out, &(state->out));
    integrator_reset(out, &(state->integ));
    arma_siso_filter_input_reset(meas, &(state->meas_deriv_filt));

    return 0;
}

int PID_update_deriv(float cmd, float meas, float meas_dot, float *out, PID_state_t *state)
{
    float cmd_filt;
    float meas_filt;
    float err_filt;

    signal_update(cmd, &cmd_filt, &(state->cmd));
    signal_update(meas, &meas_filt, &(state->cmd));
    signal_update(cmd_filt - meas_filt, &err_filt, &(state->err));

    float integ;
    integrator_update(state->Ki*err_filt, &integ, &(state->integ));
    float prop = state->Kp*err_filt;
    float ffwd = state->Kff*cmd_filt;

    float meas_dot_filt;
    signal_update(meas_dot, &meas_dot_filt, &(state->meas_dot));

    float deriv = state->Kd*meas_dot_filt;

    signal_update(ffwd + prop + integ + deriv, out, &(state->out));

    return 0;
}

int PID_update(float cmd, float meas, float *out, PID_state_t *state)
{
    float meas_dot;
    int rc = 0;
    rc = arma_siso_filter_update(meas, &meas_dot, &(state->meas_deriv_filt));
    if (rc > 0) {
        return 1;
    }

    return PID_update_deriv(cmd, meas, meas_dot, out, state);
}

