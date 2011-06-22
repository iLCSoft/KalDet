######################################################
# cmake file for building KalDet libraries
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
######################################################



# project name
PROJECT( KalDet )


# project version
SET( KalDet_VERSION_MAJOR 1 )
SET( KalDet_VERSION_MINOR 1 )
SET( KalDet_VERSION_PATCH 1 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.0 REQUIRED )
FIND_PACKAGE( KalTest 1.1.1 REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )

FOREACH( pkg Marlin MarlinUtil ROOT KalTest )

    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ENDIF()

ENDFOREACH()



### LIBRARY AND TOOLS #######################################################


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


FILE( GLOB_RECURSE ild_sources "ild/*.cc" )
FILE( GLOB_RECURSE ild_headers "ild/*.h" )
INCLUDE_DIRECTORIES( ild/common )

INSTALL( FILES ${ild_headers} DESTINATION "include/kaldet" )

ADD_SHARED_LIBRARY( KalDet ${lib_sources} ${ild_sources} )
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

