include(FetchContent)

FetchContent_Declare(
  tinygltf
  GIT_REPOSITORY https://github.com/syoyo/tinygltf.git
  GIT_TAG        v2.9.3
)

FetchContent_MakeAvailable(tinygltf)
add_library(tinygltf::tinygltf ALIAS tinygltf)

