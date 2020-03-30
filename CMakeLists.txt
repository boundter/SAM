cmake_minimum_required(VERSION 3.0)
project(sam CXX)


# version string 0.1
set(VERSION "0.1" CACHE STRING "version" FORCE)

# use c++14 standard
set(CMAKE_CXX_STANDARD 14)

# if g++ show all warnings
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall")
endif(CMAKE_COMPILER_IS_GNUCXX)

option(BUILD_TESTS "Set ON to enable tests" OFF)

# include boost
find_package(Boost 1.61.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#add_subdirectory(lib)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(sam INTERFACE)
target_include_directories(sam INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(sam ${Boost_LIBARIES})



if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
