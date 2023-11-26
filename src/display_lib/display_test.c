#include <stdio.h>
#include "pico/stdlib.h"
#include "pico/time.h"

#include "display.h"

int main() {

    stdio_init_all();

    display_init();

    uint32_t icount = 0;
    while (1) {
        
        icount = icount + 1u;
        display_page(1337, icount, 9000, 9001, true);

        sleep_ms(100);

    }

    display_close();

}
