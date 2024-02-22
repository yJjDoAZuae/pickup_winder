#include "pico/stdlib.h"

/*Encoder GPIO*/
// change these as needed
#define ENC_A	21  // Encoder phase A
#define ENC_B	20  // Encoder phase B
#define ENC_SW	22  // Encoder push botton switch

void encoder_callback(uint gpio, uint32_t event);
int get_enc_pos();
