#ifdef __arm__
    #include "pico/stdlib.h"
#else
    #include <stdlib.h>
    #include <stdbool.h>
#endif

#include "math_constants.h"

#ifndef LIMIT_H
#define LIMIT_H

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


#endif
