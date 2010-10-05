# ROOT package required
#INCLUDE( ${PROJECT_SOURCE_DIR}/FindROOT.cmake )
FIND_PACKAGE( ROOT REQUIRED )
SET( ENV{ROOTSYS} "${ROOT_HOME}" )
SET( ENV{PATH} $ENV{ROOTSYS}/bin:$ENV{PATH} )
SET( ENV{LD_LIBRARY_PATH} $ENV{ROOTSYS}/lib:$ENV{LD_LIBRARY_PATH} )
IF(APPLE)
    SET( ENV{DYLD_LIBRARY_PATH} $ENV{LD_LIBRARY_PATH} )
ENDIF(APPLE)
#-----------------------------

# ============================================================================
# prepares input headers for GEN_ROOT_DICT_SOURCES
# expects:
#   INPUT_DIR - directory to search for headers
#   HEADERS_INSTALL_DIR - if defined will intall headers to given dir
#
# returns:
#   ROOT_DICT_INPUT_HEADERS - header files found in INPUT_DIR
#       LinkDef.h should be the last header (if found) (required by rootcint)
#
MACRO( PREPARE_ROOT_DICT_HEADERS INPUT_DIR )

    FILE( GLOB ROOT_DICT_INPUT_HEADERS "${INPUT_DIR}/*.h" )

    # rootcint needs LinkDef.h as the last header !!

    #LIST( FIND ROOT_DICT_INPUT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h LINKDEF_FOUND )
    #IF( LINKDEF_FOUND )
    #    LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h )
    #    LIST( APPEND ROOT_DICT_INPUT_HEADERS LinkDef.h )
    #ENDIF( LINKDEF_FOUND )

    FILE( GLOB exclude_headers "${INPUT_DIR}/LinkDef.h" )
    IF( exclude_headers )
        LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS "${exclude_headers}" )
        IF( HEADERS_INSTALL_DIR )
            INSTALL( FILES ${ROOT_DICT_INPUT_HEADERS} DESTINATION ${HEADERS_INSTALL_DIR} )
        ENDIF()
        LIST( APPEND ROOT_DICT_INPUT_HEADERS "${exclude_headers}")
    ENDIF( exclude_headers )

    #MESSAGE( STATUS "  headers : ${ROOT_DICT_INPUT_HEADERS}  ********************  exclude headers : ${exclude_headers} " )

ENDMACRO( PREPARE_ROOT_DICT_HEADERS )

IF( NOT DEFINED ROOT_DICT_OUTPUT_DIR )
    MESSAGE( STATUS "ROOT_DICT_OUTPUT_DIR set to default: ${PROJECT_BINARY_DIR}/rootdict" )
    SET( ROOT_DICT_OUTPUT_DIR "${PROJECT_BINARY_DIR}/rootdict" )
ENDIF( NOT DEFINED ROOT_DICT_OUTPUT_DIR )

# clean generated header files with 'make clean'
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ROOT_DICT_OUTPUT_DIR}" )


IF( DEFINED ROOT_DICT_CINT_DEFINITIONS )
    MESSAGE( STATUS "ROOT_DICT_CINT_DEFINITIONS set to: ${ROOT_DICT_CINT_DEFINITIONS}" )
ENDIF( DEFINED ROOT_DICT_CINT_DEFINITIONS )

# ============================================================================
# macro for generating root dict sources with rootcint
#
# expects:
#       ROOT_DICT_INPUT_SOURCES - list of sources to generate
#       ROOT_DICT_INPUT_HEADERS - list of headers needed to generate dict sources
#       ROOT_DICT_INCLUDE_DIRS - list of include dirs to pass to rootcint -I..
#       ROOT_DICT_CINT_DEFINITIONS - extra definitions to pass to rootcint
#       ROOT_DICT_OUTPUT_DIR - where sources should be generated
#
# returns:
#       ROOT_DICT_OUTPUT_SOURCES - list containing all generated sources
# ----------------------------------------------------------------------------
MACRO( GEN_ROOT_DICT_SOURCES ROOT_DICT_INPUT_SOURCES )

#IF( NOT DEFINED ROOT_DICT_OUTPUT_DIR )
#    MESSAGE( FATAL_ERROR "ROOT_DICT_OUTPUT_DIR not defined" )
#ENDIF()

# generate dict source file
#--- need to get the include dirs in a string that rootcint understands...
set(_incs)
FOREACH(_current ${ROOT_DICT_INCLUDE_DIRS})
    SET(_incs "${_incs}\t-I${_current}")  #fg: the \t fixes a wired string expansion 
ENDFOREACH()
#MESSAGE(STATUS "--- _incs: " ${_incs} )

SET( ROOT_DICT_OUTPUT_SOURCES )
FOREACH( dict_src_filename ${ROOT_DICT_INPUT_SOURCES} )
    SET( dict_src_file ${ROOT_DICT_OUTPUT_DIR}/${dict_src_filename} )
    ADD_CUSTOM_COMMAND(
        OUTPUT  ${dict_src_file}
        COMMAND mkdir ARGS -p ${ROOT_DICT_OUTPUT_DIR}
        COMMAND rootcint
        ARGS -f "${dict_src_file}" -c ${ROOT_DICT_CINT_DEFINITIONS} ${_incs} ${ROOT_DICT_INPUT_HEADERS}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DEPENDS ${ROOT_DICT_INPUT_HEADERS}
        COMMENT "generating: ${dict_src_file}"
    )
    LIST( APPEND ROOT_DICT_OUTPUT_SOURCES ${dict_src_file} )
ENDFOREACH()

#MESSAGE( STATUS "${ROOT_DICT_OUTPUT_SOURCES}" )

ENDMACRO( GEN_ROOT_DICT_SOURCES )

# === end of macro GEN_ROOT_DICT_SOURCES =====================================

