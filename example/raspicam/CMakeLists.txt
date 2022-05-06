cmake_minimum_required(VERSION 3.1)
project(apriltag-raspicam)

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 14)

add_executable(${PROJECT_NAME} main.cpp)

include_directories("../../include")

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

# raspicam
# find_package(raspicam REQUIRED)
# include_directories(${raspicam_INCLUDE_DIRS})
# target_link_libraries(${PROJECT_NAME} ${raspicam_LIBS})
include_directories("/usr/local/include")
target_link_libraries(${PROJECT_NAME} "/usr/local/lib/libraspicam.so")
