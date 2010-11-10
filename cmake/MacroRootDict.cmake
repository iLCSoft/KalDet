IF(APPLE)
    SET( LD_LIBRARY_PATH_VAR DYLD_LIBRARY_PATH )
ELSE()
    SET( LD_LIBRARY_PATH_VAR LD_LIBRARY_PATH )
ENDIF()

IF( NOT DEFINED ROOT_DICT_OUTPUT_DIR )
    SET( ROOT_DICT_OUTPUT_DIR "${PROJECT_BINARY_DIR}/rootdict" )
ENDIF()

# clean generated header files with 'make clean'
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ROOT_DICT_OUTPUT_DIR}" )

IF( NOT ROOT_FIND_QUIETLY )
    MESSAGE( STATUS "Check for ROOT_DICT_OUTPUT_DIR: ${PROJECT_BINARY_DIR}/rootdict" )
    MESSAGE( STATUS "Check for ROOT_DICT_CINT_DEFINITIONS: ${ROOT_DICT_CINT_DEFINITIONS}" )
ENDIF()


# ============================================================================
# helper macro to prepare input headers for GEN_ROOT_DICT_SOURCES
#   sorts LinkDef.h to be the last header (required by rootcint)
#
# arguments:
#   INPUT_DIR - directory to search for headers matching *.h
#
# returns:
#   ROOT_DICT_INPUT_HEADERS - all header files found in INPUT_DIR with
#       LinkDef.h as the last header (if found)
#
# ----------------------------------------------------------------------------
MACRO( PREPARE_ROOT_DICT_HEADERS INPUT_DIR )

    FILE( GLOB ROOT_DICT_INPUT_HEADERS "${INPUT_DIR}/*.h" )
    FILE( GLOB _linkdef_hdr "${INPUT_DIR}/LinkDef.h" )

    #LIST( FIND ROOT_DICT_INPUT_HEADERS ${_linkdef_hdr} _aux )
    #IF( ${_aux} EQUAL 0 OR ${_aux} GREATER 0 )
    #    LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}" )
    #    LIST( APPEND ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}" )
    #ENDIF()

    IF( _linkdef_hdr )
        LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}" )
        LIST( APPEND ROOT_DICT_INPUT_HEADERS "${_linkdef_hdr}")
    ENDIF()

    #MESSAGE( STATUS "ROOT_DICT_INPUT_HEADERS: ${ROOT_DICT_INPUT_HEADERS}" )

ENDMACRO( PREPARE_ROOT_DICT_HEADERS )





# ============================================================================
# macro for generating root dict sources with rootcint
#
# requires following variables:
#       ROOT_DICT_INPUT_SOURCES - list of sources to generate
#       ROOT_DICT_INPUT_HEADERS - list of headers needed to generate dict sources
#           * if LinkDef.h is in the list it must be at the end !!
#       ROOT_DICT_INCLUDE_DIRS - list of include dirs to pass to rootcint -I..
#       ROOT_DICT_CINT_DEFINITIONS - extra definitions to pass to rootcint
#       ROOT_DICT_OUTPUT_DIR - where sources should be generated
#
# returns:
#       ROOT_DICT_OUTPUT_SOURCES - list containing all generated sources
# ----------------------------------------------------------------------------
MACRO( GEN_ROOT_DICT_SOURCES ROOT_DICT_INPUT_SOURCES )

    # need to prefix all include dirs with -I
    set( _dict_includes )
    FOREACH( _inc ${ROOT_DICT_INCLUDE_DIRS} )
        SET( _dict_includes "${_dict_includes}\t-I${_inc}")  #fg: the \t fixes a wired string expansion 
        #SET( _dict_includes ${_dict_includes} -I${_inc} )
    ENDFOREACH()


    SET( ROOT_DICT_OUTPUT_SOURCES )
    FOREACH( _dict_src_filename ${ROOT_DICT_INPUT_SOURCES} )
        STRING( REPLACE "/" "_" _dict_src_filename ${_dict_src_filename} )
        SET( _dict_src_file ${ROOT_DICT_OUTPUT_DIR}/${_dict_src_filename} )
        STRING( REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" _dict_hdr_file "${_dict_src_file}" )
        ADD_CUSTOM_COMMAND(
            OUTPUT  ${_dict_src_file} ${_dict_hdr_file}
            COMMAND mkdir ARGS -p ${ROOT_DICT_OUTPUT_DIR}
            COMMAND ${LD_LIBRARY_PATH_VAR}=${ROOT_LIBRARY_DIR} ${ROOT_CINT}
            ARGS -f "${_dict_src_file}" -c ${ROOT_DICT_CINT_DEFINITIONS} ${_dict_includes} ${ROOT_DICT_INPUT_HEADERS}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            DEPENDS ${ROOT_DICT_INPUT_HEADERS}
            COMMENT "generating: ${_dict_src_file} ${_dict_hdr_file}"
        )
        LIST( APPEND ROOT_DICT_OUTPUT_SOURCES ${_dict_src_file} )
    ENDFOREACH()

ENDMACRO( GEN_ROOT_DICT_SOURCES ROOT_DICT_INPUT_SOURCES )
# ============================================================================

