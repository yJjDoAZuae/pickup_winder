
aux_source_directory(. DIR_display_SRCS)

add_subdirectory(Config)
add_subdirectory(Fonts)
add_subdirectory(GUI)
add_subdirectory(LCD)

include_directories(Config)
include_directories(Fonts)
include_directories(GUI)
include_directories(LCD)

# add_library(display ${DIR_examples_SRCS})
# target_link_libraries(display PUBLIC Config Fonts GUI LCD)
