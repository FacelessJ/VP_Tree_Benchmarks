cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(vp_tree_benchmarks VERSION 0.1)

# Prevent in source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILDS ON)

include(GNUInstallDirs)

# Set up output binary directory
set(OUTPUT_DIR ${CMAKE_SOURCE_DIR}/bin)
if(UNIX)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		set(OUTPUT_DIR ${CMAKE_SOURCE_DIR}/bin/Debug_linux)
	else()
		set(OUTPUT_DIR ${CMAKE_SOURCE_DIR}/bin/Release_linux)
	endif()
endif()
message(STATUS "Output dir: ${OUTPUT_DIR}")

# Preprocessor definitions
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")
if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_WIN32 -D_CRT_SECURE_NO_WARNINGS")
elseif(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_UNIX")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})

set(SRC_DIR "src")
set(INC_DIR "include")

set(SRC_FILES
	"${SRC_DIR}/main.cpp"
	"${SRC_DIR}/basic_vp_tree.cpp"
	"${SRC_DIR}/omp_vp_tree.cpp"
	"${SRC_DIR}/async_vp_tree.cpp"
	"${SRC_DIR}/brute_force.cpp"
	"${SRC_DIR}/timing.cpp"
	)

set(INC_FILES
	"${INC_DIR}/types.hpp"
	"${INC_DIR}/basic_vp_tree.hpp"
	"${INC_DIR}/omp_vp_tree.hpp"
	"${INC_DIR}/async_vp_tree.hpp"
	"${INC_DIR}/brute_force.hpp"
	"${INC_DIR}/timing.hpp"
	)

add_executable(vp_tree_benchmarks
	${SRC_FILES}
	${INC_FILES})

source_group("src" FILES ${SRC_FILES})
source_group("include" FILES ${INC_FILES})

if(WIN32)
	if(NOT CMAKE_VERSION VERSION_LESS "3.8")
		set_target_properties(vp_tree_benchmarks PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/\$(Configuration)")
	else()
		message(SEND_ERROR "Needs CMake 3.8+")
	endif()
endif()
