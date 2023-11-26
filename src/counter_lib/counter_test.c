#include <stdio.h>
#include "pico/stdlib.h"

#include "counter.h"
#include "../motor_lib/motor.h"


void gpio_callback(uint gpio, uint32_t events) {

    if(gpio==COUNTER_GPIO_PIN && events == COUNTER_EDGE) {
        counter_gpio_edge_isr_call();
    }

}

const float speed_cmd_vec[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
int dir = 1;
const uint16_t nspeed = 10;
const uint16_t speed_mod = 40;

void main_loop(uint16_t * k)
{

    if (*k%speed_mod == 0) {
        motor_step(dir, speed_cmd_vec[*k/speed_mod], 1.0f);
    }

    counter_update();

    printf("%d count %5d, speed %5.2f, accel %5.2f\n", *k, get_count(), get_speed(), get_accel());

    (*k)++;
    if (*k > (nspeed-1)*speed_mod)
    {
        *k = 0;
        dir = -dir;
    }
}

int main()
{

    stdio_init_all();

    motor_init();

    // motor_step(1, 5.0);

    // Everything after this point happens in the PWM interrupt handler, so we
    // can twiddle our thumbs

    uint16_t k = 0;

    counter_init();
    counter_reset(0);

    gpio_set_irq_enabled_with_callback(COUNTER_GPIO_PIN, COUNTER_EDGE, true, &gpio_callback);

    while (1)
    {

        main_loop(&k);
        sleep_ms(100);
    }
}
