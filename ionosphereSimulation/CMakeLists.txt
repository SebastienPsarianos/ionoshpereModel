cmake_minimum_required(VERSION 3.16)
project(IonosphereSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SOURCES
    src/main.cpp
    src/Grid.cpp
    src/Conductance.cpp
    src/Coordinates.cpp
    src/LegacyMHDConversion.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
