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

# include Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# add_subdirectory(lib)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(sam INTERFACE)
target_include_directories(sam INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
#target_link_libraries(sam ${Boost_LIBARIES})

include(CMakePackageConfigHelpers)
set(PROJECT_CMAKE_DIR "lib/cmake/sam")
set(version_config "${PROJECT_BINARY_DIR}/sam-config-version.cmake")
set(project_config "${PROJECT_BINARY_DIR}/sam-config.cmake")
set(targets_export_name sam-targets)

write_basic_package_version_file(
  ${version_config}
  VERSION ${VERSION}
  COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/sam-config.cmake.in
  ${project_config}
  INSTALL_DESTINATION ${PROJECT_CMAKE_DIR})

export(TARGETS sam FILE
  ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

install(
  FILES ${project_config} ${version_config}
  DESTINATION ${PROJECT_CMAKE_DIR})

install(
  EXPORT ${targets_export_name}
  DESTINATION ${PROJECT_CMAKE_DIR})

install(
  TARGETS sam
  EXPORT ${targets_export_name}
  DESTINATION lib)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include
  DESTINATION ${CMAKE_INSTALL_PREFIX})


if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
