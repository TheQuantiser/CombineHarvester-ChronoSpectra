include(FindPackageHandleStandardArgs)

set(_vdt_root "${VDT_ROOT}")
if(NOT _vdt_root)
  set(_vdt_root "$ENV{CONDA_PREFIX}")
endif()

find_path(VDT_INCLUDE_DIR
  NAMES vdt/vdtMath.h
  HINTS "${_vdt_root}/include"
  NO_DEFAULT_PATH)
find_library(VDT_LIBRARY
  NAMES vdt
  HINTS "${_vdt_root}/lib"
  NO_DEFAULT_PATH)

find_package_handle_standard_args(VDT
  REQUIRED_VARS VDT_INCLUDE_DIR VDT_LIBRARY)

if(VDT_FOUND AND NOT TARGET VDT::VDT)
  add_library(VDT::VDT UNKNOWN IMPORTED)
  set_target_properties(VDT::VDT PROPERTIES
    IMPORTED_LOCATION "${VDT_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${VDT_INCLUDE_DIR}")
endif()

set(VDT_INCLUDE_DIRS "${VDT_INCLUDE_DIR}")
set(VDT_LIBRARIES "${VDT_LIBRARY}")
mark_as_advanced(VDT_INCLUDE_DIR VDT_LIBRARY)
