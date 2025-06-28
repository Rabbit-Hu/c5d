message(STATUS "Including assimp")

include(FetchContent)

FetchContent_Declare(
    assimp
    GIT_REPOSITORY https://github.com/assimp/assimp.git
    GIT_TAG v5.4.3
)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ASSIMP_INJECT_DEBUG_POSTFIX OFF CACHE BOOL "" FORCE)
set(ASSIMP_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(assimp)
# add_library(assimp::assimp ALIAS assimp)

message(STATUS "Done including Assimp")
