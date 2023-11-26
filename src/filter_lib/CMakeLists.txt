add_library(filter_lib
        filter.h
        filter.c
        )

 
add_executable(filter_test
    filter_test.c
    )

# enable serial output from the USB port
pico_enable_stdio_usb(filter_lib 1)
pico_enable_stdio_uart(filter_lib 0)

pico_enable_stdio_usb(filter_test 1)
pico_enable_stdio_uart(filter_test 0)

# pull in common dependencies
target_link_libraries(filter_lib pico_stdlib)
target_link_libraries(filter_test filter_lib pico_stdlib)

# create map/bin/hex file etc.
pico_add_extra_outputs(filter_test)

# add url via pico_set_program_url
example_auto_set_url(filter_test)
