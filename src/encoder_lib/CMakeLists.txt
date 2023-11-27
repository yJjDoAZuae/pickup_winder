add_library(encoder_lib
        encoder.h
        encoder.c
        )

add_executable(encoder_test
        encoder_test.c
        )
        
# enable serial output from the USB port
pico_enable_stdio_usb(encoder_lib 1)
pico_enable_stdio_uart(encoder_lib 0)

pico_enable_stdio_usb(encoder_test 1)
pico_enable_stdio_uart(encoder_test 0)

# pull in common dependencies
target_link_libraries(encoder_lib pico_stdlib hardware_pio)
target_link_libraries(encoder_test encoder_lib)

# create map/bin/hex file etc.
pico_add_extra_outputs(encoder_test)

# add url via pico_set_program_url
example_auto_set_url(encoder_test)
