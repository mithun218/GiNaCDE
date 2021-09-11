# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                             #
# File: cmake/FindGiNaC.cmake                                 #
#                                                             #
# - Try to find GiNaC                                         #
# Originally from Martin Velis                                #
#                                                             #
# Project name: GiNaCDE - GiNaC Differential Equation Solver  #
# Contact: Mithun Bairagi <bairagirasulpur@gmail.com>         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Version: 1.0                                                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


set( GINAC_FOUND FALSE )

find_path( GINAC_INCLUDE_DIR NAMES ginac/ginac.h)

find_library( GINAC_LIBRARY NAMES ginac)

if( GINAC_INCLUDE_DIR AND GINAC_LIBRARY )
        SET( GINAC_FOUND TRUE )
ENDIF (GINAC_INCLUDE_DIR AND GINAC_LIBRARY)

IF (GINAC_FOUND)
   IF (NOT GiNaC_FIND_QUIETLY)
      MESSAGE(STATUS "Found GiNaC: ${GINAC_LIBRARY}")
   ENDIF (NOT GiNaC_FIND_QUIETLY)
ELSE (GINAC_FOUND)
   IF (GiNaC_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GiNaC")
   ENDIF (GiNaC_FIND_REQUIRED)
ENDIF (GINAC_FOUND)

MARK_AS_ADVANCED (	GINAC_INCLUDE_DIR
					GINAC_LIBRARY
				 )
