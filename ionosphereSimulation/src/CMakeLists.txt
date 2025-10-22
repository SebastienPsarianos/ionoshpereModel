cmake_minimum_required(VERSION 3.16)
project(MyProject LANGUAGES CXX)

# Use C++17 or higher (adjust as needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add all source files
set(SOURCES
    main.cpp
    FileParser.cpp
    Conductance.cpp
)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Optionally, enable warnings
if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /permissive-)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()
