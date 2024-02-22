#include <stdio.h>
#include "pico/stdlib.h"
#include "pico/time.h"

#include "encoder.h"

int main()
{

    stdio_init_all();

	// GPIO Setup for Encoder
	gpio_init(ENC_SW);					//Initialise a GPIO for (enabled I/O and set func to GPIO_FUNC_SIO)
    gpio_set_dir(ENC_SW,GPIO_IN);
	gpio_disable_pulls(ENC_SW);

	gpio_init(ENC_A);
    gpio_set_dir(ENC_A,GPIO_IN);
	gpio_disable_pulls(ENC_A);

	gpio_init(ENC_B);
    gpio_set_dir(ENC_B,GPIO_IN);
	gpio_disable_pulls(ENC_B);

	gpio_set_irq_enabled_with_callback(ENC_SW, GPIO_IRQ_EDGE_FALL, true, &encoder_callback);
    gpio_set_irq_enabled_with_callback(ENC_A, GPIO_IRQ_EDGE_FALL, true, &encoder_callback);
	gpio_set_irq_enabled_with_callback(ENC_B, GPIO_IRQ_EDGE_FALL, true, &encoder_callback);
    gpio_set_irq_enabled_with_callback(ENC_A, GPIO_IRQ_EDGE_RISE, true, &encoder_callback);
	gpio_set_irq_enabled_with_callback(ENC_B, GPIO_IRQ_EDGE_RISE, true, &encoder_callback);

    while (1)
    {
        sleep_ms(100);
    }
}
