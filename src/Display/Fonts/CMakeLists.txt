
aux_source_directory(. DIR_Fonts_SRCS)

add_library(Fonts ${DIR_Fonts_SRCS})
target_link_libraries(Fonts PUBLIC pico_stdlib)