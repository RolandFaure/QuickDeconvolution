cmake_minimum_required(VERSION 2.8)

project(split10x)
#add_executable(${PROJECT_NAME} "main.cpp")

# Local header files here ONLY
set(TARGET_H
    tools.h
    compress.h
    build_graph.h
    create_experiments.h
    index_reads.h
    draw_10x.h
    cluster_graph.h
    generate_random_sequence.h
    measure_speed.h
    tests.h
    robin_hood.h
    read.h
    output.h
    clipp.h
   )

# Local source files here
set(TARGET_SRC
    tools.cpp
    compress.cpp
    build_graph.cpp
    create_experiments.cpp
    index_reads.cpp
    draw_10x.cpp
    cluster_graph.cpp
    generate_random_sequence.cpp
    measure_speed.cpp
    tests.cpp
    read.cpp
    output.cpp
    )

file(GLOB SOURCES "*.cpp")
add_executable(split10x ${SOURCES})

target_compile_options (split10x PRIVATE -std=c++11)
target_compile_options (split10x PRIVATE -O3)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
