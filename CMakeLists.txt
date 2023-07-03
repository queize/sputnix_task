cmake_minimum_required(VERSION 3.0)


project(sputnix_task LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)


add_executable(main ${SOURCES})
target_include_directories(main
    PUBLIC ${PROJECT_SOURCE_DIR}/include
    PUBLIC boost_1_82_0
)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/${CMAKE_BUILD_TYPE})
