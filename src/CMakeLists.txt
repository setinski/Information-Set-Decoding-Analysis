cmake_minimum_required(VERSION 3.15.0)
enable_language(CXX)
project(InformationSetDecodingAnalysis)

set(MOSEK_HEADER_DIR "C:/Program Files/Mosek/9.2/tools/platform/win64x86/h")    			 # set the path to the mosek header directory
set(MOSEK_LIB "C:/Program Files/Mosek/9.2/tools/platform/win64x86/bin/mosek64_9_2.lib") 	 # set the path to mosek library
set(FUSION_LIB "C:/Program Files/Mosek/9.2/tools/platform/win64x86/bin/fusion64_9_2.lib") 	 # set the path to fusion library
set(BOOST_DIR "C:/Program Files (x86)/boost_1_75_0")           				 				 # set the path to the boost directory

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc /std:c++latest")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

include_directories("${MOSEK_HEADER_DIR}")
include_directories("${BOOST_DIR}")

file(GLOB headers *.h)
file(GLOB sources *.c*)

add_executable("${PROJECT_NAME}" ${headers} ${sources})

set_property(TARGET "${PROJECT_NAME}" PROPERTY
				MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY 
				VS_STARTUP_PROJECT "${PROJECT_NAME}")
  
target_link_libraries("${PROJECT_NAME}" debug "${MOSEK_LIB}")
target_link_libraries("${PROJECT_NAME}" optimized "${MOSEK_LIB}")
target_link_libraries("${PROJECT_NAME}" debug "${FUSION_LIB}")
target_link_libraries("${PROJECT_NAME}" optimized "${FUSION_LIB}")
