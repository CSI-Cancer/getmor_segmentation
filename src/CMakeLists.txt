cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)

project(getmor_segmentation)

set(ITK_DIR /home/rish/code/tools/itk/itk_build)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_CXX_FLAGS "-Wall -O2")

add_executable(single_channel_mask single_channel_mask.cpp)
target_link_libraries(single_channel_mask ${ITK_LIBRARIES})

add_executable(merge_and_watershed merge_and_watershed.cpp)
target_link_libraries(merge_and_watershed ${ITK_LIBRARIES})
