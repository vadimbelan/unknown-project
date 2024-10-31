cmake_minimum_required(VERSION 3.1)

add_subdirectory(libs)

project(Reconstruction)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/opt/homebrew/opt/libomp/lib -lomp")
else()
    message(FATAL_ERROR "OpenMP not found! You need OpenMP for speedup on multicore CPUs!")
endif()

find_package(OpenCV 4.5.1 REQUIRED)

set(SOURCES
        src/reconstruction/sift/sift.cpp
        src/reconstruction/sift/sift.h
        )

set(LIBRARIES
        libutils
        ${OpenCV_LIBS}
        )

add_library(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/src)

set(TEST_SOURCES
        tests/other/other_test.cpp
        tests/other/other_test.h
        )

set(TEST_LIBRARIES
        gtest
        gtest_main
        )

add_executable(sift_main tests/sift_main.cpp ${TEST_SOURCES})
target_link_libraries(sift_main ${PROJECT_NAME} ${LIBRARIES} ${TEST_LIBRARIES})
