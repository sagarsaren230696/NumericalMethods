cmake_minimum_required(VERSION 3.10)
SET (CMAKE_CXX_COMPILER "c:/msys64/mingw64/bin/clang++.exe" CACHE STRING "C++ compiler" FORCE) # only for windows



# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
set(project_name "steeleBox")
project(${project_name})

# Get everything in the src/ directory (always), but not recursive
file(GLOB APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(${project_name} ${APP_SOURCES})

# set the directories for include files
set (APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
list (APPEND APP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
# list (APPEND APP_INCLUDE_DIRS "D:/KyawThu/Programming/KT_CPP2021/includesCPP/csv")
# list (APPEND APP_INCLUDE_DIRS "D:/KyawThu/Programming/KT_CPP2021/includesCPP/refprop")
# list (APPEND APP_INCLUDE_DIRS "D:/KyawThu/Programming/KT_CPP2021/includesCPP/numerical")

target_include_directories(${project_name} PUBLIC ${APP_INCLUDE_DIRS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (${project_name} PUBLIC Eigen3::Eigen)

# find_package(Eigen3 CONFIG REQUIRED)
# target_link_libraries(${project_name} PRIVATE Eigen3::Eigen)

# find_package(TBB REQUIRED)
# target_link_libraries(${project_name} PUBLIC TBB::tbb)
# target_include_directories(${project_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(${project_name} PUBLIC G:/Programming/KT_CPP/includesCPP/eigen3)
# target_include_directories(${project_name} PUBLIC G:/Programming/KT_CPP/includesCPP/csv)

target_compile_options(${project_name} PRIVATE -m64  -g  -lpthread  -Wall  -pedantic  -lstdc++  -std=c++20  -O3  -ltbb)

# target_compile_options(${project_name} PRIVATE -m64  -g  -lpthread  -Wall  -pedantic  -lstdc++  -std=c++20  -O3)