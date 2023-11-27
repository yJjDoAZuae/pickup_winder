#include "pico/stdlib.h"

#define DEBUG_USB 1

int display_init();

int display_page(uint32_t rev_cmd, uint32_t rev_count, 
    float speed_cmd, float speed_meas, bool direction);

int display_close();
