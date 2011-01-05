######################################################
# cmake file for building KalDet libraries
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
######################################################



# project name
PROJECT( KalDet )


# project version
SET( KALDET_VERSION_MAJOR 1 )
SET( KALDET_VERSION_MINOR 0 )
SET( KALDET_VERSION_PATCH 0 )



### SETTINGS ################################################################

# ilcsoft_default_settings.cmake needs to be in CMAKE_MODULE_PATH
INCLUDE( ilcsoft_default_settings )
#INCLUDE( build_32bit_compatible )



### DEPENDENCIES ############################################################

SET( kaldet_req_pkgs KalTest Marlin GEAR ROOT )

FOREACH( pkg ${kaldet_req_pkgs} )

    FIND_PACKAGE( ${pkg} REQUIRED )

    IF( ${pkg}_FOUND )

        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )

    ELSE()
        MESSAGE( FATAL_ERROR "required package ${pkg} not found" )
    ENDIF()

ENDFOREACH()

### LIBRARY AND TOOLS #######################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )

#-----------------------------------------------------

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

INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

#MESSAGE( STATUS "KalDet lib sources: ${lib_sources}" )

ADD_SHARED_LIBRARY( KalDet ${lib_sources} )
INSTALL_SHARED_LIBRARY( KalDet DESTINATION lib )
TARGET_LINK_LIBRARIES( KalDet ${ROOT_LIBRARIES} )

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalDet PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )


###################################################################################


# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "KalDetLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalDetLibDeps.cmake" DESTINATION lib/cmake )

