
#include "pico/stdlib.h"
#include "pico/time.h"
#include "../control_lib/control.h"
#include "pico.h"
#include "hardware/structs/timer.h"
#include "hardware/irq.h"

#define COUNTER_GPIO_PIN 4
#define COUNTER_EDGE GPIO_IRQ_EDGE_RISE

typedef struct
{
    int direction;  // 1 = CW, -1 = CCW
    uint64_t time_at_edge_us;
    // uint64_t time_last_high_us;
    uint64_t last_period_us;
    uint32_t turns_count;
    float speed_rps;
    float accel_rps2;

    arma_siso_filter_state_t log2_period_filt;
    signal_state_t speed_sig;
    arma_siso_filter_state_t accel_filt;

} counter_state_t;

int counter_init();

int counter_reset(uint32_t count);

int counter_update();

int counter_gpio_edge_isr_call();

uint32_t get_count();
float get_speed();
float get_accel();
