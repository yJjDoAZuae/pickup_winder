add_library(display_lib
        display.h
        display.c
        # ImageData.c
        # LCD_2in_test.c
        )

add_executable(display_test
        display_test.c
        )
        
# enable serial output from the USB port
pico_enable_stdio_usb(display_lib 1)
pico_enable_stdio_uart(display_lib 0)

pico_enable_stdio_usb(display_test 1)
pico_enable_stdio_uart(display_test 0)

target_include_directories(display_lib PRIVATE ../Display/Config)
target_include_directories(display_lib PRIVATE ../Display/Fonts)
target_include_directories(display_lib PRIVATE ../Display/GUI)
target_include_directories(display_lib PRIVATE ../Display/LCD)

# target_include_directories(display_test PRIVATE ../Display/Config)
# target_include_directories(display_test PRIVATE ../Display/Fonts)
# target_include_directories(display_test PRIVATE ../Display/GUI)
# target_include_directories(display_test PRIVATE ../Display/LCD)

# pull in common dependencies
target_link_libraries(display_lib pico_stdlib hardware_pwm hardware_spi hardware_i2c Config Fonts GUI LCD)
target_link_libraries(display_test display_lib)

# create map/bin/hex file etc.
pico_add_extra_outputs(display_test)

# add url via pico_set_program_url
example_auto_set_url(display_test)
