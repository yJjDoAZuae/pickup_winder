#include <stdio.h>
#include "pico/stdlib.h"
#include "pico/time.h"
//#include "hardware/irq.h"

#include "motor.h"

int main() {

    stdio_init_all();

    motor_init();

    // motor_step(1, 5.0);

    float speed_cmd_vec[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    // Everything after this point happens in the PWM interrupt handler, so we
    // can twiddle our thumbs
    
    uint16_t k = 0;
    int dir = 1;

    while (1) {
        printf("dir = %d, speed cmd = %0.1f\n", dir, speed_cmd_vec[k]);
        // motor_step(dir, speed_cmd_vec[k]);
        
        k++;
        if (k > 9) {
            k = 0;
            dir = -dir;
        }

        sleep_ms(1000);

    }
    

}
