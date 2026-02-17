cmake_minimum_required(VERSION 3.16) 
project(IonosphereSolver LANGUAGES CXX)

# --- Project Settings ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# IMPORTANT: Do NOT set CMAKE_CXX_COMPILER here. 
# Your build.sh script passes 'mpicxx' via the command line, which is the correct way.

# --- dependencies ---

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

# 3. Trilinos
# This looks for TrilinosConfig.cmake in the path provided by your build script
find_package(Trilinos REQUIRED)

# Optional: Print what we found to the terminal for debugging
message(STATUS "Trilinos Found: ${Trilinos_DIR}")
message(STATUS "Trilinos Libraries: ${Trilinos_LIBRARIES}")

# --- Main Executable ---
add_executable(${PROJECT_NAME})

# Attach sources directly to the target
target_sources(${PROJECT_NAME} PRIVATE
    src/main.cpp
    src/utils/Grid.cpp
    src/utils/Coordinates.cpp
    src/utils/LegacyMHDConversion.cpp
    src/conductance/EmpConductance.cpp
    src/conductance/utils.cpp
    src/solver/GsSolver.cpp
    src/solver/TlSolver.cpp
    src/postProcessing/EField.cpp
)

# --- Includes ---
target_include_directories(${PROJECT_NAME} PRIVATE 
    "${CMAKE_SOURCE_DIR}/include"
    ${Trilinos_INCLUDE_DIRS}      # Trilinos headers
    ${Trilinos_TPL_INCLUDE_DIRS}  # Third-party headers (MPI, BLAS, LAPACK)
)

# --- Linking ---
target_link_libraries(${PROJECT_NAME} PRIVATE 
    eigen_lib
    json_lib
    ${Trilinos_LIBRARIES}      # The main Trilinos libs
    ${Trilinos_TPL_LIBRARIES}  # The underlying libs (MPI, BLAS, etc.)
)

separate_arguments(TRILINOS_FLAGS_LIST NATIVE_COMMAND "${Trilinos_CXX_COMPILER_FLAGS}")

target_compile_options(${PROJECT_NAME} PRIVATE 
    -Wall -Wextra -Wpedantic -O3
    ${TRILINOS_FLAGS_LIST}
)
