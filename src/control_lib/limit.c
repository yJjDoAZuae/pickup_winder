#include "limit.h"

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
