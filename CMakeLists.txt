cmake_minimum_required(VERSION 3.16)
project(CAM LANGUAGES CXX)

## Check whether appropriate compiler version is used
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient gcc/g++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient clang/clang++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
endif()


## Find packages that might become necessary.
find_package(Doxygen)
find_package (Python3 COMPONENTS Interpreter Development NumPy)
find_program (CYTHON_EXECUTABLE NAMES cython cython.bat cython3)
execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${_pip_args} ${pyevtk})

## Create a target that transports all the necessary information to downstream projects
add_library(CAM INTERFACE)

## Add the correct include directories
target_include_directories(
  CAM
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

## Add the necessary C++ standard requirement
target_compile_features(
  CAM
  INTERFACE
  cxx_std_20
)

## Add the examples subdirectory
add_subdirectory(cpp_examples)


## Add a namespaced alias target for the export
add_library(CAM::CAM ALIAS CAM)

## Install targets, headers and configuration
include(GNUInstallDirs)
install(
  TARGETS CAM
  EXPORT CAM-config
)

install(
  EXPORT CAM-config
  NAMESPACE CAM::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CAM
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

## Write a configuration file for Cython
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake_cython.cfg.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_cython.cfg
)

## Build doxygen if found.
if(DOXYGEN_FOUND)
  set(DOXYGEN_TAB_SIZE 2)
  set(DOXYGEN_TOC_INCLUDE_HEADINGS 0)
  set(DOXYGEN_EXTRACT_ALL YES)
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_EXTRACT_STATIC YES)
  set(DPXYGEN_STRICT_PROTO_MATCHING YES)
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_MACRO_EXPANSION YES)
  set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
  set(DOXYGEN_CALL_GRAPH YES)
  set(DOXYGEN_DOT_MULTI_TARGETS NO)
  set(DOXYGEN_EXAMPLE_PATH ${CMAKE_SOURCE_DIR}/examples)
  set(DOXYGEN_FILTER_PATTERNS *.m=../submodules/doxymatlab.git/m2cpp.pl)
  set(DOXYGEN_EXCLUDE_PATTERNS */submodules/*)
  doxygen_add_docs(
    doxygen
    ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT Building Doxygen documentation...
  )
endif()


## Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
