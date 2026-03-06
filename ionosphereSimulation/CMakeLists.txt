cmake_minimum_required(VERSION 3.16) 
project(IonosphereSolver LANGUAGES CXX)

# --- Project Settings ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


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
find_package(Trilinos REQUIRED)

# 4. Build Exec
add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE
    src/main.cpp
    src/utils/Grid.cpp
    src/utils/Coordinates.cpp
    src/utils/LegacyMHDConversion.cpp
    src/conductance/EmpConductance.cpp
    src/conductance/utils.cpp
    src/conductance/SolarModel.cpp
    src/conductance/DipoleModel.cpp
    src/solver/GsSolver.cpp
    src/solver/TlSolver.cpp
    src/postProcessing/EField.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    "${CMAKE_SOURCE_DIR}/include"
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    eigen_lib
    json_lib
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
)

separate_arguments(TRILINOS_FLAGS_LIST NATIVE_COMMAND "${Trilinos_CXX_COMPILER_FLAGS}")

target_compile_options(${PROJECT_NAME} PRIVATE 
    -Wall -Wextra -Wpedantic -O3
    ${TRILINOS_FLAGS_LIST}
)
