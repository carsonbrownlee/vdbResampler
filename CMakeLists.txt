cmake_minimum_required (VERSION 2.6)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project (VDBResampler)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
find_package (OpenVDB REQUIRED)
find_package (TBB REQUIRED COMPONENTS tbb)
find_package (IlmBase REQUIRED COMPONENTS Half)
message( "openvdb dirs: " ${OpenVDB_INCLUDE_DIRS})
include_directories(${OpenVDB_INCLUDE_DIRS})
add_executable(vdbResampler vdbResampler.cpp)
message( "openvdb libs: " ${OpenVDB_LIBRARIES})
target_link_libraries(vdbResampler 
   ${OpenVDB_LIBRARIES} ${Tbb_LIBRARIES} ${IlmBase_LIBRARIES})
