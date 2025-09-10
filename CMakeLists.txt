cmake_minimum_required(VERSION 3.10)
project(pacejka)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

include(FetchContent)
FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)
FetchContent_MakeAvailable(matplotplusplus)

add_library(pacejka pacejka.cpp)

add_executable(run_example run_example.cpp)
target_link_libraries(run_example pacejka matplot)