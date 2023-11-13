cmake_minimum_required(VERSION 3.10.0)
project(dpm VERSION 0.1.0 LANGUAGES C CXX)

# Set C++ Standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

include(CTest)
enable_testing()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Automatically include all source files
file(GLOB SOURCES "src/*.cpp")
# Print all source files
message(STATUS "SOURCES: ${SOURCES}")

add_executable(dpm main.cpp ${SOURCES})

# Link against HDF5 libraries
target_link_libraries(dpm ${HDF5_CXX_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
