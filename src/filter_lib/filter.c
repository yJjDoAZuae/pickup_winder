#include "math.h"
#include <stdio.h>

#include "filter.h"

float clip(float in, float vmin, float vmax)
{
    return fmin(fmax(in, vmin), vmax);
}

int arma_buffer_init(arma_buffer_t *buff)
{
    for (int k = 0; k < ARMA_SISO_MAX_ORDER+1; k++) {
        buff->k[k] = 0.0;
    }
    return 0;
}

int arma_siso_filter_init(arma_siso_filter_state_t *state)
{
    const float tol = 1e-6;

    state->n = 0;

    arma_buffer_init(&(state->u));
    arma_buffer_init(&(state->y));
    arma_buffer_init(&(state->a));
    arma_buffer_init(&(state->b));

    // we ignore the first entry in den and set it to 1 in all cases
    state->a.k[0] = 1.0f;
    
    // pass through with DC gain of 1.0
    state->b.k[0] = 1.0f;

    return 0;
}


int tustin_1(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz)
{
    const float tol = 1.0e-6;

    float K = 2.0f/dt;
    float coeff_denom = den->k[0]*K + den->k[1];

    // printf("tustin_1: den->k[0] = %0.3f, den->k[1] = %0.3f, K = %0.3f\n", den->k[0], den->k[1], K);

    if (fabs(coeff_denom) < tol) {
        return 1;
    }
    
    numz->k[0] = (num->k[0]*K + num->k[1])/coeff_denom;
    numz->k[1] = (-num->k[0]*K + num->k[1])/coeff_denom;

    denz->k[0] = 1.0f;
    denz->k[1] = (-den->k[0]*K + den->k[1])/coeff_denom;

    return 0;
}

int tustin_2(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz)
{
    const float tol = 1.0e-6;

    float K = 2.0f/dt;
    float coeff_denom = den->k[0]*K*K + den->k[1]*K + den->k[2];

    if (fabs(coeff_denom) < tol) {
        return 1;
    }
    
    numz->k[0] = (     num->k[0]*K*K + num->k[1]*K +     num->k[2])/coeff_denom;
    numz->k[1] = (-2.0*num->k[0]*K*K               + 2.0*num->k[2])/coeff_denom;
    numz->k[2] = (     num->k[0]*K*K - num->k[1]*K +     num->k[2])/coeff_denom;

    denz->k[0] = 1.0f;
    denz->k[1] = (-2.0*den->k[0]*K*K               + 2.0*den->k[2])/coeff_denom;
    denz->k[2] = (     den->k[0]*K*K - den->k[1]*K +     den->k[2])/coeff_denom;

    return 0;
}

int arma_siso_lowpass_1(float dt, float tau, arma_siso_filter_state_t * state)
{

    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // HACK: zero order hold realization
    // ydot = -1/tau y + 1/tau u
    // yk+1 - yk = (-1/tau yk + 1/tau u) * dt
    // yk+1 = (1 - dt/tau) xk + dt/tau u

    // numz.k[0] = dt/tau;
    // numz.k[1] = 0.0f;
    // denz.k[0] = 1.0f;
    // denz.k[1] = -(1-dt/tau);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 0.0f;
    num.k[1] = 1/tau;
    den.k[0] = 1.0f;
    den.k[1] = 1/tau;

    // printf("num[0] = %0.3f, num[1] = %0.3f, den[0] = %0.3f, den[1] = %0.3f\n", num.k[0], num.k[1], den.k[0], den.k[1]);

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        printf("tustin returned %d\n", rc);
        return rc;
    }

    // printf("numz[0] = %0.3f, numz[1] = %0.3f, denz[0] = %0.3f, denz[1] = %0.3f\n", numz.k[0], numz.k[1], denz.k[0], denz.k[1]);

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_lowpass_2(float dt, float wn, float zeta, float tau_zero, arma_siso_filter_state_t *state)
{

    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 0.0f;
    num.k[1] = tau_zero;
    num.k[2] = wn*wn;
    den.k[0] = 1.0f;
    den.k[1] = 2.0*zeta*wn;
    den.k[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_notch_2(float dt, float wn, float zeta_den, 
    float zeta_num, arma_siso_filter_state_t *state)
{

    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 1.0f;
    num.k[1] = 2.0*zeta_num*wn;
    num.k[2] = wn*wn;
    den.k[0] = 1.0f;
    den.k[1] = 2.0*zeta_den*wn;
    den.k[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_highpass_1(float dt, float tau, arma_siso_filter_state_t * state)
{
    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 1/tau;
    num.k[1] = 0.0f;
    den.k[0] = 1.0f;
    den.k[1] = 1/tau;

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_highpass_2(float dt, float wn, float zeta, float c_zero, arma_siso_filter_state_t *state)
{
    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 1.0f;
    num.k[1] = c_zero*2.0*zeta*wn;
    num.k[2] = 0.0f;
    den.k[0] = 1.0f;
    den.k[1] = 2.0*zeta*wn;
    den.k[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }
    
    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_deriv_1(float dt, float tau, arma_siso_filter_state_t * state)
{
    arma_buffer_t numz;
    arma_buffer_t denz;

    arma_buffer_init(&numz);
    arma_buffer_init(&denz);

    // Tustin (bilinear) realization
    arma_buffer_t num;
    arma_buffer_t den;

    arma_buffer_init(&num);
    arma_buffer_init(&den);

    num.k[0] = 1/tau;
    num.k[1] = 0.0f;
    den.k[0] = 1.0f;
    den.k[1] = 1/tau;

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_filter_coef_update(arma_buffer_t *num,
                                 arma_buffer_t *den, uint8_t n, arma_siso_filter_state_t *state)
{

    arma_buffer_init(&(state->a));
    arma_buffer_init(&(state->b));

    if (n > ARMA_SISO_MAX_ORDER)
    {
        return 1;
    }

    state->n = n;

    for (int k = 0; k < state->n + 1; k++)
    {
        state->b.k[k] = num->k[k];
        state->a.k[k] = den->k[k];
    }

    // we ignore the first entry in den and set it to 1 in all cases
    state->a.k[0] = 1.0f;

    return 0;
}

int arma_siso_filter_input_reset(float in, arma_siso_filter_state_t * state)
{
    const float tol = 1e-6;

    arma_buffer_init(&(state->u));
    arma_buffer_init(&(state->y));

    // TRICKY: potential access to uninitialized data
    // for zeroth order (static gain) case.
    float numsum = 0.0f;
    float densum = 0.0f;

    for (int k = 0; k < state->n+1; k++) {
        numsum += state->b.k[k];
        densum += state->a.k[k];
    }

    // We need an arbitrary finite DC gain to accomodate 
    // band pass, high pass, and lead/lag filters filters
    // TODO: Filter stability test?

    float dcgain;
    // Test for infinite DC gain.
    // Pure integrators should not be implemented
    // using an ARMA filter.
    if (fabs(densum) > tol) {
        dcgain = numsum/densum;

        printf("dcgain = %0.3f\n", dcgain);
    } else {
        return 2;
    }

    for (int k = 0; k < ARMA_SISO_MAX_ORDER; k++) {
        state->u.k[k] = in;
        state->y.k[k] = in*dcgain;
    }

    return 0;
}

int arma_siso_filter_update(float in, float *out, arma_siso_filter_state_t *state)
{

    *out = state->b.k[0]*in;

    // NOTE: implicit state->a.k[0] == 1 because
    // we're assigning state->y.k[0] without a coefficient

    // update the output
    for (int k = 1; k < state->n+1; k++) {
        // TRICKY: because we haven't rolled the buffers yet
        // the input and output buffer indices are one less
        // than the coefficient indices
        *out += state->b.k[k] * state->u.k[k-1];
        *out -= state->a.k[k] * state->y.k[k-1];
    }

    // roll the buffers.  Always keep the whole buffer in case we change filter coefficients.
    for (int k = ARMA_SISO_MAX_ORDER; k > 0; k--) {
        state->y.k[k] = state->y.k[k-1];
        state->u.k[k] = state->u.k[k-1];
    }
    state->y.k[0] = *out;
    state->u.k[0] = in;

    return 0;
}

int limit_init(limit_state_t *state)
{

    state->dt = 1.0f;

    state->min = 0.0f;
    state->max = 0.0f;

    state->lower = false;
    state->upper = false;
    state->y = 0.0f;

    return 0;
}

int val_limit_update(float in, float *out, limit_state_t *state)
{

    if (state->max < state->min) {
        return 1;
    }

    state->y = in;

    state->lower = state->y < state->min;
    state->upper = state->y > state->max;

    if (state->lower) {
        state->y = state->min;
    }

    if (state->upper) {
        state->y = state->max;
    }

    *out = state->y;

    return 0;
}

int rate_limit_reset(float in, limit_state_t *state)
{

    state->lower = false;
    state->upper = false;
    state->y = in;

    return 0;
}

int rate_limit_update(float in, float * out, limit_state_t *state)
{

    if (state->max < state->min) {
        return 1;
    }

    if (state->dt <= 0.0f) {
        return 2;
    }

    float del = in - state->y;

    state->lower = del < state->min/state->dt;
    state->upper = del > state->max/state->dt;

    if (state->lower) {
        del = state->min/state->dt;
    }

    if (state->upper) {
        del = state->max/state->dt;
    }

    state->y += del;

    *out = state->y;

    return 0;
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

    arma_siso_deriv_1(dt, state->deriv_tau, &(state->meas_deriv_filt));

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

