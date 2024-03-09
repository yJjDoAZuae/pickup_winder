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

int arma_buffer_init(arma_buffer_t *buff)
{
    memset(*buff, 0, (ARMA_SISO_MAX_ORDER+1)*sizeof(float));
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
    state->a[0] = 1.0f;
    
    // pass through with DC gain of 1.0
    state->b[0] = 1.0f;

    return 0;
}


int tustin_1(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz)
{
    const float tol = 1.0e-6;

    float K = 2.0f/dt;
    float coeff_denom = *den[0]*K + *den[1];

    // printf("tustin_1: den->k[0] = %0.3f, den->k[1] = %0.3f, K = %0.3f\n", den->k[0], den->k[1], K);

    if (fabs(coeff_denom) < tol) {
        return 1;
    }
    
    *numz[0] = (*num[0]*K + *num[1])/coeff_denom;
    *numz[1] = (-*num[0]*K + *num[1])/coeff_denom;

    *denz[0] = 1.0f;
    *denz[1] = (-*den[0]*K + *den[1])/coeff_denom;

    return 0;
}

int tustin_2(arma_buffer_t * num, arma_buffer_t * den, float dt, arma_buffer_t * numz, arma_buffer_t * denz)
{
    const float tol = 1.0e-6;

    float K = 2.0f/dt;
    float coeff_denom = *den[0]*K*K + *den[1]*K + *den[2];

    if (fabs(coeff_denom) < tol) {
        return 1;
    }
    
    *numz[0] = (     *num[0]*K*K + *num[1]*K +     *num[2])/coeff_denom;
    *numz[1] = (-2.0*(*num)[0]*K*K               + 2.0*(*num)[2])/coeff_denom;
    *numz[2] = (     *num[0]*K*K - *num[1]*K +     *num[2])/coeff_denom;

    *denz[0] = 1.0f;
    *denz[1] = (-2.0*(*den)[0]*K*K               + 2.0*(*den)[2])/coeff_denom;
    *denz[2] = (     *den[0]*K*K - *den[1]*K +     *den[2])/coeff_denom;

    return 0;
}

int arma_siso_lowpass_1_design(float dt, float tau, arma_siso_filter_state_t * state)
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

    num[0] = 0.0f;
    num[1] = 1/tau;
    den[0] = 1.0f;
    den[1] = 1/tau;

    // printf("num[0] = %0.3f, num[1] = %0.3f, den[0] = %0.3f, den[1] = %0.3f\n", num.k[0], num.k[1], den.k[0], den.k[1]);

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        // printf("tustin returned %d\n", rc);
        return rc;
    }

    // printf("numz[0] = %0.3f, numz[1] = %0.3f, denz[0] = %0.3f, denz[1] = %0.3f\n", numz.k[0], numz.k[1], denz.k[0], denz.k[1]);

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_lowpass_2_design(float dt, float wn, float zeta, float tau_zero, arma_siso_filter_state_t *state)
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

    num[0] = 0.0f;
    num[1] = tau_zero;
    num[2] = wn*wn;
    den[0] = 1.0f;
    den[1] = 2.0*zeta*wn;
    den[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_notch_2_design(float dt, float wn, float zeta_den, 
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

    num[0] = 1.0f;
    num[1] = 2.0*zeta_num*wn;
    num[2] = wn*wn;
    den[0] = 1.0f;
    den[1] = 2.0*zeta_den*wn;
    den[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_highpass_1_design(float dt, float tau, arma_siso_filter_state_t * state)
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

    num[0] = 1/tau;
    num[1] = 0.0f;
    den[0] = 1.0f;
    den[1] = 1/tau;

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_highpass_2_design(float dt, float wn, float zeta, float c_zero, arma_siso_filter_state_t *state)
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

    num[0] = 1.0f;
    num[1] = c_zero*2.0*zeta*wn;
    num[2] = 0.0f;
    den[0] = 1.0f;
    den[1] = 2.0*zeta*wn;
    den[2] = wn*wn;

    int rc = tustin_2(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }
    
    return arma_siso_filter_coef_update(&numz, &denz, 2, state);
}

int arma_siso_deriv_1_design(float dt, float tau, arma_siso_filter_state_t * state)
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

    num[0] = 1/tau;
    num[1] = 0.0f;
    den[0] = 1.0f;
    den[1] = 1/tau;

    int rc = tustin_1(&num, &den, dt, &numz, &denz);
    if (rc!=0) {
        return rc;
    }

    return arma_siso_filter_coef_update(&numz, &denz, 1, state);
}

int arma_siso_filter_coef_update(arma_buffer_t *num,
                                 arma_buffer_t *den, int n, arma_siso_filter_state_t *state)
{

    arma_buffer_init(&(state->a));
    arma_buffer_init(&(state->b));

    if (n < 0 || n > ARMA_SISO_MAX_ORDER)
    {
        return 1;
    }

    state->n = n;

    for (int k = 0; k < state->n + 1; k++)
    {
        state->b[k] = *num[k];
        state->a[k] = *den[k];
    }

    // we ignore the first entry in den and set it to 1 in all cases
    state->a[0] = 1.0f;

    return 0;
}

int arma_siso_filter_input_reset(float in, arma_siso_filter_state_t * state)
{
    const float tol = 1e-6;

    arma_buffer_init(&(state->u));
    arma_buffer_init(&(state->y));

    float dcgain;
    int rc = arma_siso_filter_dcgain(&dcgain, state);

    if (rc == 0) {
        for (int k = 0; k < ARMA_SISO_MAX_ORDER; k++) {
            state->u[k] = in;
            state->y[k] = in*dcgain;
        }
    }

    return 0;
}

int arma_siso_filter_output_reset(float out, arma_siso_filter_state_t * state)
{
    const float tol = 1e-6;

    arma_buffer_init(&(state->u));
    arma_buffer_init(&(state->y));

    float dcgain;
    int rc = arma_siso_filter_dcgain(&dcgain, state);

    if (rc == 0 && fabs(dcgain) > tol) {
        for (int k = 0; k < ARMA_SISO_MAX_ORDER; k++) {
            state->y[k] = out;
            state->u[k] = out/dcgain;
        }
    } else {
        return 1;
    }

    return 0;
}

int arma_siso_filter_dcgain(float * dcgain, arma_siso_filter_state_t *state)
{
    const float tol = 1e-6;

    // TRICKY: potential access to uninitialized data
    // for zeroth order (static gain) case.
    float numsum = 0.0f;
    float densum = 0.0f;

    for (int k = 0; k < state->n+1; k++) {
        numsum += state->b[k];
        densum += state->a[k];
    }

    // We need an arbitrary finite DC gain to accomodate 
    // band pass, high pass, and lead/lag filters filters
    // TODO: Filter stability test?
    *dcgain = 0.0f;

    // Test for infinite DC gain.
    // Pure integrators should not be implemented
    // using an ARMA filter.
    if (fabs(densum) > tol) {
        *dcgain = numsum/densum;

        // printf("dcgain = %0.3f\n", dcgain);
    } else {
        return 1;
    }

    return 0;
}

int arma_siso_filter_update(float in, float *out, arma_siso_filter_state_t *state)
{

    *out = state->b[0]*in;

    // NOTE: implicit state->a.k[0] == 1 because
    // we're assigning state->y.k[0] without a coefficient

    // update the output
    for (int k = 1; k < state->n+1; k++) {
        // TRICKY: because we haven't rolled the buffers yet
        // the input and output buffer indices are one less
        // than the coefficient indices
        *out += state->b[k] * state->u[k-1];
        *out -= state->a[k] * state->y[k-1];
    }

    // roll the buffers.  Always keep the whole buffer in case we change filter coefficients.
    for (int k = ARMA_SISO_MAX_ORDER; k > 0; k--) {
        state->y[k] = state->y[k-1];
        state->u[k] = state->u[k-1];
    }
    state->y[0] = *out;
    state->u[0] = in;

    return 0;
}

// Init a FIR buffer to zeros
int fir_buffer_init(fir_buffer_t *buff)
{
    memset(*buff, 0, FIR_SISO_MAX_ORDER*sizeof(fir_value_t));
    return 0;
}

// Initialize a FIR filter
int fir_siso_filter_init(fir_siso_filter_state_t *state)
{

    int rc = 0;

    state->n = 0;

    rc = fir_buffer_init(&state->b);
    if (rc!=0) {
        return rc;
    }

    rc = fir_buffer_init(&state->u);
    if (rc!=0) {
        return rc;
    }
    
    return rc;
}


// Reset a filter based on a constant input value
int fir_siso_filter_dcgain(fir_value_t * dcgain, fir_siso_filter_state_t *state)
{
    long long tmp_out = 0;

    if (state->n > FIR_SISO_MAX_ORDER){
        return 1;
    }

    for (int k = 0; k < state->n; k++) {
        tmp_out += state->b[k];
    }

    *dcgain = tmp_out/state->n;

    return 0;
}

// Reset a filter based on a constant input value
int fir_siso_filter_input_reset(fir_value_t in, fir_siso_filter_state_t *state)
{
    if (state->n > FIR_SISO_MAX_ORDER){
        return 1;
    }

    for (int k = state->n-1; k >= 0; k--) {
        state->u[k] = in;
    }

    return 0;
}

// Reset a filter based on a constant output value
int fir_siso_filter_output_reset(fir_value_t out, fir_siso_filter_state_t *state)
{
    long long sum_b = 0;

    if (state->n > FIR_SISO_MAX_ORDER){
        return 1;
    }

    for (int k = 0; k < state->n; k++) {
        sum_b += state->b[k];
    }

    fir_value_t in = 0;
    if (sum_b != 0) {
        in = out*state->n / sum_b;
    } else {
        in = 0;
    }

    for (int k = 0; k < state->n; k++) {
        state->u[k] = in;
    }

    if (sum_b == 0) {
        return 2;
    }

    return 0;
}

// Step a filter forward one iteration
int fir_siso_filter_update(fir_value_t in, fir_value_t *out, fir_siso_filter_state_t *state)
{
    long long tmp_out = 0;

    if (state->n > FIR_SISO_MAX_ORDER){
        return 1;
    }

    for (int k = state->n-1; k > 0; k--) {
        state->u[k] = state->u[k-1];
        tmp_out += state->u[k] * state->b[k];
    }
    state->u[0] = in;
    tmp_out += state->u[0]*state->b[0];

    *out = tmp_out/state->n;

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

