######################################################
# cmake file for building KalDet libraries
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
######################################################



# project name
PROJECT( KalDet )


# project version
SET( KalDet_VERSION_MAJOR 1 )
SET( KalDet_VERSION_MINOR 8 )
SET( KalDet_VERSION_PATCH 0 )

cmake_policy(SET CMP0008 NEW)


### DEPENDENCIES ###########################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES streamlog REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

 
#FG: switch back to depend on Marlin (still needed by lctpc/lp1jgem/EXTPCMeasLayer.cxx)
FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( GEAR REQUIRED )
FIND_PACKAGE( KalTest 1.2 REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )

#NB: need ROOT before KalTest (needed on case sensitive systems and classes that differ only in capitalization)
FOREACH( pkg Marlin ROOT KalTest LCIO GEAR )

  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
  ENDIF()

ENDFOREACH()

### LIBRARY AND TOOLS ######################################################
 

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )


# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( ${ROOT_DICT_MACROS_FILE} )



# build KalDet library
SET( lib_input_dirs gen kern lctpc/lp1jgem )

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${KalTest_INCLUDE_DIRS} )

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )



FOREACH( lib_input_dir ${lib_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )

    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )

    INSTALL_DIRECTORY( ${lib_input_dir}/ DESTINATION "include/kaldet"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE
    )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )

    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

#MESSAGE( "KalDet lib sources: ${lib_sources}" )


INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )


# ild directory is treated separately as it does not require any root dictionaries

FILE( GLOB_RECURSE ild_sources "ild/*.cc" )
FILE( GLOB_RECURSE lib_headers "ild/*.h" )

INCLUDE_DIRECTORIES( ild/common )

INSTALL( FILES ${lib_headers} DESTINATION "include/kaldet" )

# needed for adding header files to xcode project
IF(CMAKE_GENERATOR MATCHES "Xcode")

  FOREACH( lib_input_dir ${lib_input_dirs} ) 

    FILE( GLOB_RECURSE lib_input_headers "${lib_input_dir}/*.h" )
    LIST( APPEND  lib_headers ${lib_input_headers} ) 

  ENDFOREACH()   

  ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${lib_sources} ${ild_sources} ${lib_headers} )

ELSE()
  ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${lib_sources} ${ild_sources} )
ENDIF()      


INSTALL_SHARED_LIBRARY( KalDet DESTINATION lib )

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalDet PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )



# display some variables and write them to cache
DISPLAY_STD_VARIABLES()


# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( KalDetConfig.cmake KalDetConfigVersion.cmake KalDetLibDeps.cmake )

