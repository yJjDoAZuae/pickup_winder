add_library(control_lib
        control.h
        math_constants.h
        control.c
        )

 
add_executable(control_test
    control_test.c
    )

# enable serial output from the USB port
pico_enable_stdio_usb(control_lib 1)
pico_enable_stdio_uart(control_lib 0)

pico_enable_stdio_usb(control_test 1)
pico_enable_stdio_uart(control_test 0)

# pull in common dependencies
target_link_libraries(control_lib pico_stdlib)
target_link_libraries(control_test control_lib pico_stdlib)

# create map/bin/hex file etc.
pico_add_extra_outputs(control_test)

# add url via pico_set_program_url
example_auto_set_url(control_test)
