#include <stdio.h>
#include "math.h"
#include "pico/time.h"
#include "hardware/irq.h"
#include "hardware/gpio.h"

#include "counter.h"

static counter_state_t state;

const float small_speed = 0.5f;

int counter_init() {

    state.direction = 1;
    state.time_at_edge_us = 0;
    // state.time_last_high_us = 0;
    state.last_period_us = 0;
    state.turns_count = 0;
    state.speed_rps = 0.0f;
    state.accel_rps2 = 0.0f;

    const float dt = 0.1;
    const float period_tau = 0.5;
    const float deriv_tau = 1.0;

    int rc = 0;

    rc = arma_siso_filter_init(&(state.log2_period_filt));
    if (rc > 0) {
        return 1;
    }

    rc = signal_init(&(state.speed_sig));
    if (rc > 0) {
        return 2;
    }

    rc = arma_siso_filter_init(&(state.accel_filt));
    if (rc > 0) {
        return 3;
    }

    gpio_init(COUNTER_GPIO_PIN);
    gpio_set_dir(COUNTER_GPIO_PIN, false);
    gpio_pull_up(COUNTER_GPIO_PIN);

    arma_siso_lowpass_1(dt, period_tau, &(state.log2_period_filt));
    signal_init(&(state.speed_sig));

    state.speed_sig.vlim.min = 0.0f;
    state.speed_sig.vlim.max = 15.0f;
    state.speed_sig.rlim.min = -5.0f;
    state.speed_sig.rlim.max = 5.0f;
    state.speed_sig.rlim.dt = dt;

    arma_siso_deriv_1(dt, deriv_tau, &(state.accel_filt));

    return 0;
};

// TODO: should this also reset the speed output?
int counter_reset(uint32_t count) {
    state.turns_count = count;
    state.time_at_edge_us = time_us_64();
    // state.time_last_high_us = time_us_64();
    state.last_period_us = 0;

    arma_siso_filter_input_reset(0.0f, &state.log2_period_filt);
};

int counter_update() {

    uint64_t this_ticks_us = time_us_64();

    float inst_log2_period_us;

    state.last_period_us = (uint64_t)clip(state.last_period_us, 1e6/small_speed, 1e6/15.0);

    // if the period is increasing and we're already past the expected edge, use the latest duration
    if (this_ticks_us < state.time_at_edge_us + state.last_period_us) {
        inst_log2_period_us = log2f(fmax(state.last_period_us, 1.0f));
    } else {
        inst_log2_period_us = log2f(fmax(this_ticks_us - state.time_at_edge_us, 1.0f));
    }

    float log2_period_us;

    // filter log2 of last_period_us to get a filtered period
    arma_siso_filter_update(inst_log2_period_us, &log2_period_us, &(state.log2_period_filt));

    float speed_inst_rps = 1000000.0f/(powf(2.0f, log2_period_us));

    // calculate speed
    signal_update(speed_inst_rps, &(state.speed_rps), &(state.speed_sig));

    if (state.speed_rps < small_speed)
    {
        state.speed_rps = 0.0f;
    }

    // calculate acceleration
    arma_siso_filter_update(state.speed_rps, &(state.accel_rps2), &(state.accel_filt));

};

int counter_gpio_edge_isr_call() {

    uint64_t this_ticks_us;
    
    this_ticks_us = time_us_64();

    if (this_ticks_us > state.time_at_edge_us + state.last_period_us/5) {

        state.last_period_us = this_ticks_us - state.time_at_edge_us;
        state.turns_count += 1;
        state.time_at_edge_us = this_ticks_us;

    }

};

float get_speed() {
    return state.speed_rps;
}

float get_accel() {
    return state.accel_rps2;
}

uint32_t get_count() {
    return state.turns_count;
}
