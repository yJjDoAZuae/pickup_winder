add_executable(counter_test
        counter_test.c
        )

add_library(counter_lib
        counter.h
        counter.c
        )

target_include_directories(counter_test PRIVATE ../motor_lib)

# enable serial output from the USB port
pico_enable_stdio_usb(counter_lib 1)
pico_enable_stdio_uart(counter_lib 0)

pico_enable_stdio_usb(counter_test 1)
pico_enable_stdio_uart(counter_test 0)

# pull in common dependencies
target_link_libraries(counter_lib pico_stdlib hardware_pwm hardware_spi hardware_i2c motor_lib)
target_link_libraries(counter_test counter_lib)

# create map/bin/hex file etc.
pico_add_extra_outputs(counter_test)

# add url via pico_set_program_url
example_auto_set_url(counter_test)
