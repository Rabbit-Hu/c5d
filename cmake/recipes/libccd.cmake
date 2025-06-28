# Use FetchContent to download and build libccd

include(FetchContent)

FetchContent_Declare(
  libccd
  GIT_REPOSITORY https://github.com/danfis/libccd.git
  GIT_TAG        v2.1
)

set(ENABLE_DOUBLE_PRECISION ON)

FetchContent_MakeAvailable(libccd)

# target_compile_options(ccd PUBLIC -DENABLE_DOUBLE_PRECISION=ON)

add_library(libccd::libccd ALIAS ccd)


