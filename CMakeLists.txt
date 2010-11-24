############################################
# cmake file for building KalDet libraries
# @author Jan Engels, DESY
############################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( KalDet )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )

# project version
SET( KALDET_MAJOR_VERSION 1 )
SET( KALDET_MINOR_VERSION 0 )
SET( KALDET_PATCH_LEVEL 0 )

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

##---------------- 32/64 bit issues ---------------------------------------
#OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
##INCLUDE( CheckTypeSize )
##CHECK_TYPE_SIZE( "void*" SIZE_OF_VOID_P )
#
#IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
#    MESSAGE( STATUS "32 bit architecture detected" )
#ENDIF()
#
#IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#    MESSAGE( STATUS "64 bit architecture detected" )
#
#    IF( BUILD_32BIT_COMPATIBLE )
#        IF( COMMAND SET_PROPERTY )
#            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
#            SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32" )
#        ELSE()
#            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
#        ENDIF()
#        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
#    ELSE()
#        MESSAGE( STATUS "Build in native 64 bit mode" )
#    ENDIF()
#ENDIF()
##-------------------------------------------------------------------------

# library *nix style versioning
SET( KALDET_SOVERSION "${KALDET_MAJOR_VERSION}.${KALDET_MINOR_VERSION}" )
SET( KALDET_VERSION "${KALDET_SOVERSION}.${KALDET_PATCH_LEVEL}" )

# cmake output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

###################################################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )

#-----------------------------------------------------

# macros for generating root dict sources with rootcint
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )

###################################################################################


SET( kaldet_req_pkgs KalTest Marlin GEAR ROOT )

FOREACH( pkg ${kaldet_req_pkgs} )

    IF( ${pkg}_HOME AND NOT ${pkg}_DIR )
        SET( ${pkg}_DIR ${${pkg}_HOME} )
    ENDIF()

    FIND_PACKAGE( ${pkg} REQUIRED )

    IF( ${pkg}_FOUND )

        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ELSE()
        MESSAGE( FATAL_ERROR "required package ${pkg} not found" )
    ENDIF()

ENDFOREACH()

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

    INSTALL( DIRECTORY ${lib_input_dir}/ DESTINATION "include/kaldet"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )

    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

#MESSAGE( STATUS "KalDet lib sources: ${lib_sources}" )

ADD_LIBRARY( KalDet ${lib_sources} )

ADD_CUSTOM_TARGET( lib DEPENDS KalDet )

TARGET_LINK_LIBRARIES( KalDet ${ROOT_LIBRARIES} )

SET_TARGET_PROPERTIES( KalDet PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${KALDET_VERSION}
    SOVERSION ${KALDET_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
)

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalDet PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )

# install library
INSTALL( TARGETS KalDet DESTINATION lib PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)


###################################################################################


# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall
    "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
)

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/KalDetConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/KalDetConfig.cmake" @ONLY
)

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalDetConfig.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
#IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
#ENDIF()
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
#SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
#    "Set to OFF to build static libraries" FORCE )
#SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
#    "Where to install KALDET" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "KalDetBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalDetBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "KalDetLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalDetLibDeps.cmake" DESTINATION lib/cmake )

