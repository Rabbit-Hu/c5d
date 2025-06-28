add_executable(datagen)

add_subdirectory(distances)

target_link_libraries(
    datagen
    PRIVATE
    c5d
)

set(DATAGEN_SOURCES
    DataPoint.hpp
    DataConfig.hpp
    Dataset.hpp
    datagen.cpp
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Source Files" FILES ${DATAGEN_SOURCES})
target_sources(datagen PRIVATE ${DATAGEN_SOURCES})

add_executable(compare)

target_link_libraries(
    compare
    PRIVATE
    c5d
)

set(COMPARE_SOURCES
    compare.cpp
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Source Files" FILES ${COMPARE_SOURCES})
target_sources(compare PRIVATE ${COMPARE_SOURCES})
