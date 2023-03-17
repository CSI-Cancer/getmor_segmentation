cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)

project(getmor_segmentation)

set(ITK_DIR /home/rish/code/tools/itk/itk_build)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_CXX_FLAGS "-Wall -O2")

add_executable(ck_segment ck_segment.cpp)
target_link_libraries(ck_segment ${ITK_LIBRARIES})

add_executable(single_channel_mask single_channel_mask.cpp)
target_link_libraries(single_channel_mask ${ITK_LIBRARIES})

add_executable(merge_and_watershed merge_and_watershed.cpp)
target_link_libraries(merge_and_watershed ${ITK_LIBRARIES})

add_executable(feature_extraction feature_extraction.cpp)
target_link_libraries(feature_extraction ${ITK_LIBRARIES})

add_executable(filter_features filter_features.cpp)
target_link_libraries(filter_features ${ITK_LIBRARIES})

add_executable(cell_image cell_image.cpp)
target_link_libraries(cell_image ${ITK_LIBRARIES})

add_executable(gain gain.cpp)
target_link_libraries(gain ${ITK_LIBRARIES})
