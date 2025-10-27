cmake_minimum_required(VERSION 3.16)
project(IonosphereSolver LANGUAGES CXX)

# Use C++17 or higher (adjust as needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add all source files
set(SOURCES
    src/main.cpp
    src/FileParser.cpp
    src/Conductance.cpp
)


add_executable(${PROJECT_NAME} ${SOURCES})



target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
