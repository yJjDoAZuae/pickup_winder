add_library(motor_lib
        motor.h
        motor.c
        )

add_executable(motor_test
        motor_test.c
        )

# enable serial output from the USB port
pico_enable_stdio_usb(motor_test 1)
pico_enable_stdio_uart(motor_test 0)
    
# pull in common dependencies
target_link_libraries(motor_lib pico_stdlib hardware_pwm control_lib counter_lib)
target_link_libraries(motor_test motor_lib)

# create map/bin/hex file etc.
pico_add_extra_outputs(motor_test)

# add url via pico_set_program_url
example_auto_set_url(motor_test)
