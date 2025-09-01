project(nt VERSION 1.1.0 LANGUAGES CXX)

add_library(nt INTERFACE)

target_include_directories(nt
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(nt INTERFACE cxx_std_20)

add_library (nt::nt ALIAS nt)