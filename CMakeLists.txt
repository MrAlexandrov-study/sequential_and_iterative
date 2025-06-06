cmake_minimum_required(VERSION 3.14)
project(project_template CXX)

set(CMAKE_CXX_STANDARD 20)

# for coverage analyzing
# if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
# endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Boost COMPONENTS range algorithm)

file(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/main.cpp)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# add_library(Project_templateLib INTERFACE)
# target_include_directories(Project_templateLib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/main.cpp)
# target_link_libraries(${PROJECT_NAME} Project_templateLib)

enable_testing()
add_subdirectory(${CMAKE_SOURCE_DIR}/tests)

set(CTEST_BINARY_DIRECTORY ${CMAKE_BINARY_DIR})
