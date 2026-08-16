cmake_minimum_required(VERSION 3.20)
project(streaming_simplification
        VERSION 0.1.0
        DESCRIPTION "Streaming simplification of trajectory streams"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SSK_BUILD_TESTS "Build unit tests" ON)
option(SSK_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

# src/algo is the include root, so headers are included by their path within the
# library: #include "io/trajectory.hpp".
add_library(ssk
  src/algo/io/json.cpp
  src/algo/io/trajectory.cpp
  src/algo/simplify/dots.cpp
  src/algo/simplify/squish.cpp
)
target_include_directories(ssk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/algo)

if(MSVC)
  target_compile_options(ssk PRIVATE /W4)
else()
  target_compile_options(ssk PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion)
  # The geometry built on top of this is specified over IEEE doubles with explicit
  # tolerances; fast-math style reassociation would silently invalidate its predicates.
  target_compile_options(ssk PRIVATE -ffp-contract=off)
endif()
if(SSK_WARNINGS_AS_ERRORS)
  target_compile_options(ssk PRIVATE $<IF:$<BOOL:${MSVC}>,/WX,-Werror>)
endif()

if(SSK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
