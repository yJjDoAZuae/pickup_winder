#include <stdio.h>
#include "hardware/pio.h"
#include "hardware/gpio.h"

#include "encoder.h"


static int enc_pos = 0;  // Encoder position state

// bool getBit(uint32_t field, uint8_t bitNum) {
//     return (field & (0x1 << (bitNum - 1)));
// }

// CW rotation 
//             ________        
//        ___R1        F1________  Phase A
//                  ________    
//        ________R2        F2___  Phase B	
// state:  0    1    2    3    0
// wave:  0x0  0x1  0x3  0x2  0x0
// NOTE: Phase A leads Phase B for CW rotation

// CCW Rotation
//                  ________
//        ________R2        F2___  Phase A	
//             ________
//        ___R1        F1________  Phase B
// state:  0    3    2    1    0
// wave:  0x0  0x2  0x3  0x1  0x0
// NOTE: Phase B leads Phase A for CCW rotation

// TODO: for future generalization
#define WLEN 4
static const uint8_t waveform[WLEN] = {0, 1, 3, 2};

void encoder_callback(uint gpio, uint32_t event)
{
    
    uint32_t gpio_state = 0;

    // TODO: generalize this to more than 2 encoder inputs and arbitrary waveforms


    // TODO: do we still need to use this bitmask
    // or is it better to just test each gpio separately?
    // gpio_state = (gpio_get_all() >> 10) & 0b0111;  	// get all GPIO them mask out all but bits 10, 11, 12
                                                    // This will need to change to match which GPIO pins are being used.

    bool encA_state = gpio_get(ENC_A);
    bool encB_state = gpio_get(ENC_B);
    bool encSW_state = gpio_get(ENC_SW);
    
    static bool ccw_fall = false;  // bool used when falling edge that implies CCW rotation is detected
    static bool cw_fall = false;   // bool used when falling edge that implies CW rotation is detected
    static bool ccw_rise = false;  // bool used when falling edge that implies CCW rotation is detected
    static bool cw_rise = false;   // bool used when falling edge that implies CW rotation is detected
    
    // uint8_t enc_value = 0;
    // enc_value = (gpio_state & 0x03);

    if (gpio == ENC_A && event == GPIO_IRQ_EDGE_FALL) 
    {
        if (!cw_fall && encB_state)
            // detect initial edge of CW rotation
            cw_fall = true;

        if (ccw_fall && !encB_state)
        {
            // CCW rotation was already detected with an initial B fall, second fall on A confirms it
            // increment the position counter
            cw_fall = false;
            ccw_fall = false;
            cw_rise = false;
            ccw_rise = false;

            enc_pos -= 1;
            printf("CCW pos: %d\n", enc_pos);
        }

    }	

    if (gpio == ENC_B && event == GPIO_IRQ_EDGE_FALL) 
    {
        if (!ccw_fall && encA_state)
            // detect initial edge of CCW rotation
            ccw_fall = true;

        if (cw_fall && !encA_state)
        {
            // CCW rotation was already detected with an initial B fall, second fall on A confirms it
            // increment the position counter
            cw_fall = false;
            ccw_fall = false;
            cw_rise = false;
            ccw_rise = false;

            enc_pos += 1;
            printf("CW  pos: %d\n", enc_pos);
        }

    }
    
    if (gpio == ENC_A && event == GPIO_IRQ_EDGE_RISE) 
    {
        if (!cw_rise && !encB_state)
            // detect initial edge of CW rotation
            cw_rise = true; 

        if (ccw_rise && encB_state)
        {
            // CCW rotation was already detected with an initial B rise, second rise on A confirms it
            // increment the position counter
            cw_fall = false;
            ccw_fall = false;
            cw_rise = false;
            ccw_rise = false;

            enc_pos -= 1;
            printf("CCW pos: %d\n", enc_pos);
        }

    }	

    if (gpio == ENC_B && event == GPIO_IRQ_EDGE_RISE) 
    {
        if (!ccw_rise && !encA_state)
            // detect initial edge of CCW rotation
            ccw_rise = true;

        if (cw_rise && encA_state)
        {
            // CW rotation was already detected with an initial A rise, second rise on B confirms it
            // increment the position counter
            cw_fall = false;
            ccw_fall = false;
            cw_rise = false;
            ccw_rise = false;

            enc_pos += 1;
            printf("CW  pos: %d\n", enc_pos);
        }

    }

    if (gpio == ENC_SW && event == GPIO_IRQ_EDGE_FALL) 
    {
        printf("Boop!\n");
    }	
}

int get_enc_pos()
{
    return enc_pos;
}
