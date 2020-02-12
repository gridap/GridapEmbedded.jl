#################################################################
# HEADER
#################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
PROJECT(qhullwrapper C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)
SET(LIB ${PROJECT_NAME})

#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# DEFINE PATHS
#################################################################

SET(CMAKE_PATH ${CMAKE_SOURCE_DIR}/CMake)
SET(SRC_PATH ${CMAKE_SOURCE_DIR}/src)
SET(CMAKE_MODULE_PATH "${CMAKE_PATH}/Modules")
SET(LIB_PATH ${SRC_PATH}/lib)
SET(EXAMPLES_PATH ${SRC_PATH}/examples)
SET(TESTS_PATH ${SRC_PATH}/tests)

#################################################################
# BUILD PATHS
#################################################################

SET(ROOT_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${ROOT_OUTPUT_PATH}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${ROOT_OUTPUT_PATH}/bin)

#################################################################
# CONFIGURATION TYPES & BUILD MODE & BUILD_TESTS
#################################################################

SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT ${PROJECT_NAME}_ENABLE_TESTS)
    OPTION(${PROJECT_NAME}_ENABLE_TESTS "Enable/disable tests compilation" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_TESTS)

IF(NOT ${PROJECT_NAME}_ENABLE_EXAMPLES)
    OPTION(${PROJECT_NAME}_ENABLE_EXAMPLES "Enable/disable examples compilation" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_EXAMPLES)

#################################################################
# Add external projects
#################################################################

UNSET(${PROJECT_NAME}_EXTERNAL_PROJECTS CACHE)
SET(EXTERNAL_PROJECTS QHULL) # List of Fempar external projects

# Loop over all (required) external projects
# ${EXT_PROJ}_DIR is used as a clue to find every external project
FOREACH(EXT_PROJ ${EXTERNAL_PROJECTS})
    IF(NOT ${EXT_PROJ}_DIR AND ${PROJECT_NAME}_THIRDPARTY_DIR)
        # If the directory of the external project was not defined, build it from ThirdParty directory 
        SET(${EXT_PROJ}_DIR ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_THIRDPARTY_DIR}/${EXT_PROJ})
        IF(NOT EXISTS ${${EXT_PROJ}_DIR} OR NOT IS_DIRECTORY ${${EXT_PROJ}_DIR})
            GET_FILENAME_COMPONENT(${EXT_PROJ}_DIR ${${PROJECT_NAME}_THIRDPARTY_DIR}/${EXT_PROJ} ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR})
        ENDIF()
    ENDIF() 

    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_DIR} CACHE STRING "${EXT_PROJ} root directory")
    FIND_PACKAGE(${EXT_PROJ} REQUIRED)                    # Look for the external project
    SET(${PROJECT_NAME}_EXTERNAL_PROJECTS ${${PROJECT_NAME}_EXTERNAL_PROJECTS} ${${EXT_PROJ}_LIBS}) 
    MESSAGE("LIBS: ${${PROJECT_NAME}_EXTERNAL_PROJECTS}")   
    INCLUDE_DIRECTORIES(${${EXT_PROJ}_INCLUDES})          # Include directory of the external project
    SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${${EXT_PROJ}_INCLUDES})
ENDFOREACH ()


#################################################################
# FLAGS depend on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)

SET(OPENMP_FLAGS "")
IF(DEFINED ${PROJECT_NAME}_ENABLE_OPENMP AND ${PROJECT_NAME}_ENABLE_OPENMP)
    FIND_PACKAGE( OpenMP )
    IF(${OPENMP_FOUND})
        SET(OPENMP_FLAGS ${OpenMP_Fortran_FLAGS})
    ENDIF()
ENDIF()

message(STATUS "COMPILER INFO: ${CMAKE_C_COMPILER_ID} - ${C_COMPILER_NAME}")

IF (${CMAKE_C_COMPILER_ID} STREQUAL "GNU" OR C_COMPILER_NAME MATCHES "gcc*")
    # gcc 
#    set (CMAKE_C_FLAGS " " CACHE STRING "" FORCE)
ELSEIF (${CMAKE_C_COMPILER_ID} STREQUAL "Intel" OR C_COMPILER_NAME MATCHES "icc*")
    # icc
#    set (CMAKE_C_FLAGS " " CACHE STRING "" FORCE)
ENDIF ()

message (STATUS "CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
message (STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS})
message (STATUS "CMAKE_C_FLAGS_RELEASE: " ${CMAKE_C_FLAGS_RELEASE})
message (STATUS "CMAKE_C_FLAGS_DEBUG: " ${CMAKE_C_FLAGS_DEBUG})

#################################################################
# ENABLE TESTING
#################################################################

ENABLE_TESTING()
INCLUDE(CTest)

#################################################################
# DYNAMIC LIBRARY
#################################################################
SET(BUILD_SHARED_LIBS True)

#################################################################
# ADD SOURCE SUBDIRS
#################################################################

ADD_SUBDIRECTORY(${LIB_PATH})
IF(${PROJECT_NAME}_ENABLE_TESTS)
  ADD_SUBDIRECTORY(${TESTS_PATH})
ENDIF()

IF(${PROJECT_NAME}_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(${EXAMPLES_PATH})
ENDIF()

#################################################################
# CREATE CONFIG FILE
#################################################################

CONFIGURE_FILE(${CMAKE_PATH}/${PROJECT_NAME}Config.cmake.in 
               ${ROOT_OUTPUT_PATH}/${PROJECT_NAME}Config.cmake
               @ONLY)

