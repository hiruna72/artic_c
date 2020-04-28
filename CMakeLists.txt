# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.")
endif()

# choose the architecture 
set(DEPLOY_PLATFORM "unknown" CACHE STRING "choose the instruction set architecture")
set_property(CACHE DEPLOY_PLATFORM PROPERTY STRINGS unknown x86 arm64-v8a armeabi-v7a)
IF (${DEPLOY_PLATFORM} STREQUAL "unknown")
	message(FATAL_ERROR "choose the DEPLOY_PLATFORM") 
	return() # This is to stop proceeding further and to stop opencv getting set to the default ANDROID_ABI
ENDIF()


# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 OLD)

# Specify the project, and set the use of c++ and c
project(Artic_c LANGUAGES CXX C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fPIC -rdynamic -Wall -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -rdynamic -Wall -O2 -std=c++11 ")


include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/multiIntersectBed)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty)

set(trim src/trim.cpp)
set(main src/main.cpp)
set(mask src/mask.cpp)
set(multiintersct_main src/multiIntersectBed/multiIntersectBedMain.cpp)
set(multiintersct src/multiIntersectBed/multiIntersectBed.cpp)
set(bedfile src/multiIntersectBed/bedFile/bedFile.cpp)
set(genomefile src/multiIntersectBed/GenomeFile/GenomeFile.cpp)
set(gzstream src/multiIntersectBed/gzstream/gzstream.C)
set(filetype src/multiIntersectBed/fileType/fileType.cpp)
# Build a shared lib
add_executable(artic_c ${main} ${trim} ${mask} ${multiintersct_main} ${multiintersct} ${bedfile} ${genomefile} ${gzstream} ${filetype})
#add_library(artic_c SHARED ${common} ${trim})

# htslib dependency
set (htslib_shared "${PROJECT_SOURCE_DIR}/3rdparty/${DEPLOY_PLATFORM}/libhts.so")

IF (${DEPLOY_PLATFORM} STREQUAL "x86")
	target_link_libraries(artic_c -lz -ldl -lm -lpthread ${htslib_shared})
ELSE()
	target_link_libraries(artic_c -lz -ldl -lm ${htslib_shared})
ENDIF()
