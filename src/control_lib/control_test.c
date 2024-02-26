#include <stdio.h>
#include <math.h>
#include "pico/stdlib.h"
#include "pico/time.h"

#include "control.h"

int main() {

    stdio_init_all();

    float u = 0.0;
    float y = 0.0;
    float dt = 0.1;
    float tau = 2.0;
    float wn = 2.0f;
    float zeta = 0.707;

    arma_siso_filter_state_t state;

    // tustin realization of a first order discrete lowpass IIR filter
    // arma_siso_lowpass_1_init(u, dt, tau, &state);
    arma_siso_filter_init(&state);
    arma_siso_lowpass_2_design(dt, wn, zeta, 0.0f, &state);
    arma_siso_filter_input_reset(0.0f, &state);

    u = 1.0;

    float y2 = 0.0f;
    float y2_z1 = 0.0f;
    float y2dot = 0.0f;
    float u_z1 = 0.0f;

    uint32_t icount = 0;
    while (1) {
        
        if (icount > (int)ceil(10/dt)) {
            icount = 0;
            arma_siso_filter_input_reset(0.0f, &state);
            y2 = 0.0f;
        }

        printf("b[0] = %0.3f, b[1] = %0.3f, a[0] = %0.3f, a[1] = %0.3f\n", state.b.k[0], state.b.k[1], state.a.k[0], state.a.k[1]);
        printf("u[0] = %0.3f, u[1] = %0.3f, y[0] = %0.3f, y[1] = %0.3f\n", state.u.k[0], state.u.k[1], state.y.k[0], state.y.k[1]);

        arma_siso_filter_update(u, &y, &state);

        // 1st order ZOH low pass
        // y2 = (1 - dt/tau)*y2 + dt/tau*u;

        // 2nd order ZOH low pass
        float tmp = y2 + y2dot*dt;
        y2dot = y2dot + (-2.0f*zeta*wn*y2dot -wn*wn*y2 + wn*wn*u)*dt;
        
        // y2_z1 = y2;
        y2 = tmp;
        // u_z1 = u;

        printf("icount = %d, u = %0.2f, y = %0.3f, y2 = %0.3f, y-y2 = %0.3f\n", icount, u, y, y2, y-y2);

        icount = icount + 1u;

        sleep_ms(100);

    }

}
