# -----------------------------------------------------------------------------
#
# Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.19.3)
project(angiogenesis)

# BioDynaMo curretly uses the C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Use BioDynaMo in this project.
find_package(BioDynaMo REQUIRED)

# Add include directory for unibn octree
include_directories("${CMAKE_SOURCE_DIR}/third_party/unibn-octree")

# See UseBioDynaMo.cmake in your BioDynaMo build folder for details.
# Note that BioDynaMo provides gtest header/libraries in its include/lib dir.
include(${BDM_USE_FILE})

# Add the shared library libXMLIO.so in ROOT_LIBRARY_DIR to the required libraries
set(BDM_REQUIRED_LIBRARIES ${BDM_REQUIRED_LIBRARIES} "${ROOT_LIBRARY_DIR}/libXMLIO.so")

# Export compile commands for Sonar Lint
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Consider all files in src/ for BioDynaMo simulation.
include_directories("src")
file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cc)

bdm_add_executable(${CMAKE_PROJECT_NAME}
                   HEADERS ${PROJECT_HEADERS}
                   SOURCES ${PROJECT_SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

# Consider all files in test/ for GoogleTests.
include_directories("test")
file(GLOB_RECURSE TEST_SOURCES test/*.cc)
file(GLOB_RECURSE TEST_HEADERS test/*.h)

bdm_add_test(${CMAKE_PROJECT_NAME}-test
             SOURCES ${TEST_SOURCES}
             HEADERS ${TEST_HEADERS}
             LIBRARIES ${BDM_REQUIRED_LIBRARIES} ${CMAKE_PROJECT_NAME})

# Copy the folder test/data to the build directory
file(COPY test/data DESTINATION ${CMAKE_BINARY_DIR})
