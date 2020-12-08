# - Find nauty
# Find the native NAUTY headers and libraries.
#
#  NAUTY_INCLUDE_DIRS - where to find nauty.h, etc.
#  NAUTY_LIBRARIES    - List of libraries when using nauty.
#  NAUTY_FOUND        - True if nauty found.
#

find_path(NAUTY_INCLUDE_DIR nauty/nauty.h)

find_library(NAUTY_LIBRARY NAMES nauty)

set(NAUTY_LIBRARIES ${NAUTY_LIBRARY} )
set(NAUTY_INCLUDE_DIRS ${NAUTY_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NAUTY_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NAUTY DEFAULT_MSG NAUTY_LIBRARY NAUTY_INCLUDE_DIR )

mark_as_advanced(NAUTY_INCLUDE_DIR NAUTY_LIBRARY )
