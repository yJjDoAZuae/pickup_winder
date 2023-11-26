
#include "pico/stdlib.h"

void motor_init(void);

void motor_step(int direction, float speed, float dt);

uint16_t getPWM(float speed);
