cmake_minimum_required(VERSION 3.15.0)
enable_language(CXX)

#if(MSVC OR MSYS OR MINGW)
#    #
#endif()

#if(APPLE)
#    # for MacOS X or iOS, watchOS, tvOS (since 3.10.3)
#endif()

#if(UNIX AND NOT APPLE)
#    # for Linux, BSD, Solaris, Minix
#endif()

project(InformationSetDecodingAnalysis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(MOSEK_HEADER_DIR "C:/Program Files/Mosek/9.2/tools/platform/win64x86/h")    			 # set the path to the mosek header directory
set(MOSEK_LIB "C:/Program Files/Mosek/9.2/tools/platform/win64x86/bin/mosek64_9_2.lib") 	 # set the path to mosek library
set(FUSION_LIB "C:/Program Files/Mosek/9.2/tools/platform/win64x86/bin/fusion64_9_2.lib") 	 # set the path to fusion library
set(BOOST_DIR "C:/Program Files (x86)/boost_1_75_0")           				 				 # set the path to the boost directory

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

if(MSVC)
    add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
    add_compile_options("-Wall" "-Wextra" "-Werror" "$<$<CONFIG:RELEASE>:-O3>")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options("-stdlib=libc++")
    else()
        # nothing special for gcc at the moment
    endif()
endif()

include_directories("${MOSEK_HEADER_DIR}")
include_directories("${BOOST_DIR}")

file(GLOB headers *.h)
file(GLOB sources *.c*)

add_executable("${PROJECT_NAME}" ${headers} ${sources})

if(MSVC OR MSYS OR MINGW)
    set_property(TARGET "${PROJECT_NAME}" PROPERTY
				MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY 
				VS_STARTUP_PROJECT "${PROJECT_NAME}")
endif()

target_link_libraries("${PROJECT_NAME}" "${MOSEK_LIB}")
target_link_libraries("${PROJECT_NAME}" "${FUSION_LIB}")

