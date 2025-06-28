set(SOURCES
    distances.hpp
    ee_classify.hpp
    ee_classify.cpp
    vt_classify.hpp
    vt_classify.cpp
    vv_distance.hpp
    vv_distance.cpp
    ve_distance.hpp
    ve_distance.cpp
    vt_distance.hpp
    vt_distance.cpp
    ee_distance.hpp
    ee_distance.cpp
    distance_utils.hpp
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Source Files" FILES ${SOURCES})
target_sources(c5d PRIVATE ${SOURCES})