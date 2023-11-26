#include "math.h"
#include "pico/stdlib.h"
#include "hardware/pwm.h"
#include "pico/time.h"
#include "hardware/irq.h"
#include "motor.h"
#include "../filter_lib/filter.h"

const uint PIN_A = 2U;
const uint PIN_B = 3U;

#define NUMDAT 23U

#define MIN_SPEED 0.0f
#define MAX_SPEED 15.0f

const float pwmVec[NUMDAT]  = { 0,  8000, 10000, 12000, 15000, 18000, 
                21000, 24000, 27000, 30000, 33000, 36000, 
                39000, 42000, 45000, 48000, 51000, 54000, 
                57000, 60000, 63000, 65000, 65535 };
    
const float speedVec[NUMDAT] = { 0.0,  0.01,  0.97,  1.49,  2.24,  2.94, 
                 3.65,  4.33,   5.0,  5.65,  6.30,  6.95, 
                 7.64,  8.28,  8.92,  9.54, 10.20, 10.98, 
                11.62, 12.40, 13.15, 13.65, 13.79 };

const uint16_t pwmMin = 0U;
const uint16_t pwmMax = 65535U;

static int direction;
static float speed_cmd;
static float speed_meas;
static uint16_t motor_cmd;

static PID_state_t speed_loop;


void on_pwm_wrap() 
{

    // Clear the interrupt flag that brought us here
    pwm_clear_irq(pwm_gpio_to_slice_num(PIN_A));

    // TODO: do we have to do something with PIN_B's irq slice?  Same slice?
    // TODO: should we be filtering the motor command here?

    if (direction>0) {
        pwm_set_gpio_level(PIN_A, motor_cmd);
        pwm_set_gpio_level(PIN_B, pwmMin);
    } else {
        pwm_set_gpio_level(PIN_A, pwmMin);
        pwm_set_gpio_level(PIN_B, motor_cmd);
    }
 }

uint16_t getPWM(float in_speed)
{
    
    in_speed = fmax(fmin(in_speed, speedVec[NUMDAT-1]), speedVec[0]);
    
    uint16_t ks = 0;
    while (ks<NUMDAT-2 && in_speed > speedVec[ks+1]) {
        ks++;
    }

    float fCmd = pwmVec[ks] + (pwmVec[ks+1] - pwmVec[ks])/(speedVec[ks+1] - speedVec[ks])*(in_speed - speedVec[ks]);
    fCmd = fmin(fmax(fCmd, pwmMin), pwmMax);
    uint16_t cmd = (int16_t)fCmd;

    return cmd;
}

void motor_init(void) 
{

    direction = 1;
    speed_cmd = 5.0f;

    PID_init(&speed_loop);

    float speed_out = 0.0f;
    speed_meas = 0.0f;
    PID_reset(speed_cmd, speed_meas, speed_out, 0.0f, &speed_loop);
    motor_cmd = getPWM(speed_out);

    // gpio_init(PIN2); // TODO: not needed?
    // gpio_set_dir(PIN2, GPIO_OUT); TODO: not needed?

    gpio_set_function(PIN_A, GPIO_FUNC_PWM);
    gpio_set_function(PIN_B, GPIO_FUNC_PWM);

    // Figure out which slice we just connected to the LED pin
    uint slice_1 = pwm_gpio_to_slice_num(PIN_A);

    pwm_set_gpio_level(PIN_B, pwmMin);

    // Mask our slice's IRQ output into the PWM block's single interrupt line,
    // and register our interrupt handler
    pwm_clear_irq(slice_1);
    pwm_set_irq_enabled(slice_1, true);
    irq_set_exclusive_handler(PWM_IRQ_WRAP, on_pwm_wrap);
    irq_set_enabled(PWM_IRQ_WRAP, true);
    
   // Get some sensible defaults for the slice configuration. By default, the
    // counter is allowed to wrap over its maximum range (0 to 2**16-1)
    pwm_config config = pwm_get_default_config();
    
    // Set divider, reduces counter clock to sysclock/this value
    pwm_config_set_clkdiv(&config, 4.f);

    // Load the configuration into our PWM slice, and set it running.
    pwm_init(slice_1, &config, true);

}

void motor_step(int in_direction, float in_speed_cmd, float dt)
{
    direction = in_direction;

    // TODO: should we be filtering/rate limiting the speed command here?

    PID_set_dt(dt, &speed_loop);

    speed_loop.cmd.vlim.min = 0.0f;
    speed_loop.cmd.vlim.max = 15.0f;
    speed_loop.cmd.rlim.min = -0.25f;
    speed_loop.cmd.rlim.max = 0.25f;

    speed_loop.meas.vlim.min = 0.0f;
    speed_loop.meas.vlim.max = 15.0f;
    speed_loop.meas.rlim.min = -3.0f;
    speed_loop.meas.rlim.max = 3.0f;

    arma_siso_lowpass_1(dt, 1.0f, &(speed_loop.meas.filter));

    speed_loop.err.vlim.min = -5.0f;
    speed_loop.err.vlim.max = 5.0f;
    speed_loop.err.rlim.min = -10.0f;
    speed_loop.err.rlim.max = 10.0f;

    speed_loop.integ.vlim.min = -10.0f;
    speed_loop.integ.vlim.max = 10.0f;

    speed_loop.integ.dt = dt;
    speed_loop.integ.aw.k[0] = speed_meas < MIN_SPEED;
    speed_loop.integ.aw.k[1] = speed_meas > MAX_SPEED;
    speed_loop.integ.aw.dir[0] = AW_LOWER;
    speed_loop.integ.aw.dir[1] = AW_UPPER;

    speed_loop.out.vlim.min = 0.0f;
    speed_loop.out.vlim.max = 15.0f;
    speed_loop.out.rlim.min = -0.25f;
    speed_loop.out.rlim.max = 0.25f;

    arma_siso_lowpass_1(dt, 1.0f, &(speed_loop.out.filter));

    speed_loop.Kff = 1.0f;
    speed_loop.Kp = 0.1f;
    speed_loop.Ki = speed_loop.Kp/2.0f;
    speed_loop.Kd = 0.1f;
    

    speed_cmd = in_speed_cmd;

    float speed_out = speed_cmd;
    // PID_update(speed_cmd, speed_meas, &speed_out, &speed_loop);

    motor_cmd = getPWM(speed_out);
}
