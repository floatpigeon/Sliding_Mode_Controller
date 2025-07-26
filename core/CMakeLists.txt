cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

project(SlidingModeController VERSION 1.0.0 LANGUAGES CXX)

include_directories(./include)

file(GLOB_RECURSE sources CONFIGURE_DEPENDS 
src/*.cpp 
)

add_executable(${PROJECT_NAME}_example ${sources})