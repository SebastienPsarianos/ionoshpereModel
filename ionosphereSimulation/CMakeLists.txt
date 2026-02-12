cmake_minimum_required(VERSION 3.16)
project(IonosphereSolver LANGUAGES CXX)

# --- Project Settings ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 1. Eigen
add_library(eigen_lib INTERFACE)
target_include_directories(eigen_lib INTERFACE 
    "${CMAKE_SOURCE_DIR}/external/eigen-5.0.0"
)

# 2. JSON
add_library(json_lib INTERFACE)
target_include_directories(json_lib INTERFACE 
    "${CMAKE_SOURCE_DIR}/external/nlohmann"
)

# --- Main Executable ---
add_executable(${PROJECT_NAME})

# Attach sources directly to the target
target_sources(${PROJECT_NAME} PRIVATE
    src/main.cpp
    src/utils/Grid.cpp
    src/utils/Coordinates.cpp
    src/utils/LegacyMHDConversion.cpp
    src/conductance/Conductance.cpp
    src/conductance/utils.cpp
    src/solver/Solver.cpp
    src/postProcessing/EField.cpp
)

# --- Include Directories & Linking ---

# 1. Internal Includes (Your own headers)
target_include_directories(${PROJECT_NAME} PRIVATE 
    "${CMAKE_SOURCE_DIR}/include"
)

# 2. Link External Dependencies
# This automatically pulls in the include paths defined above for Eigen and JSON.
target_link_libraries(${PROJECT_NAME} PRIVATE 
    eigen_lib
    json_lib
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -O3)
