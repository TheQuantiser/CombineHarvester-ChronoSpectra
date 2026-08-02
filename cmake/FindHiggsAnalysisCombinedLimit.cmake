include(FindPackageHandleStandardArgs)

set(_cl_root "${HiggsAnalysisCombinedLimit_ROOT}")
if(NOT _cl_root)
  set(_cl_root "$ENV{CONDA_PREFIX}")
endif()
if(NOT _cl_root)
  set(_cl_root "")
endif()

find_path(HiggsAnalysisCombinedLimit_INCLUDE_DIR
  NAMES HiggsAnalysis/CombinedLimit/interface/Combine.h
  HINTS "${_cl_root}/include"
  NO_DEFAULT_PATH)
find_library(HiggsAnalysisCombinedLimit_LIBRARY
  NAMES HiggsAnalysisCombinedLimit
  HINTS "${_cl_root}/lib"
  NO_DEFAULT_PATH)
find_file(HiggsAnalysisCombinedLimit_PCM
  NAMES libHiggsAnalysisCombinedLimit_rdict.pcm
  HINTS "${_cl_root}/lib"
  NO_DEFAULT_PATH)
find_file(HiggsAnalysisCombinedLimit_ROOTMAP
  NAMES libHiggsAnalysisCombinedLimit.rootmap
  HINTS "${_cl_root}/lib"
  NO_DEFAULT_PATH)

find_package_handle_standard_args(HiggsAnalysisCombinedLimit
  REQUIRED_VARS
    HiggsAnalysisCombinedLimit_INCLUDE_DIR
    HiggsAnalysisCombinedLimit_LIBRARY
    HiggsAnalysisCombinedLimit_PCM
    HiggsAnalysisCombinedLimit_ROOTMAP)

if(HiggsAnalysisCombinedLimit_FOUND)
  get_filename_component(HiggsAnalysisCombinedLimit_PREFIX
    "${HiggsAnalysisCombinedLimit_INCLUDE_DIR}/.." ABSOLUTE)
  if(NOT TARGET HiggsAnalysisCombinedLimit::HiggsAnalysisCombinedLimit)
    add_library(HiggsAnalysisCombinedLimit::HiggsAnalysisCombinedLimit UNKNOWN IMPORTED)
    set_target_properties(HiggsAnalysisCombinedLimit::HiggsAnalysisCombinedLimit
      PROPERTIES
        IMPORTED_LOCATION "${HiggsAnalysisCombinedLimit_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${HiggsAnalysisCombinedLimit_PREFIX}/include")
  endif()
  set(HiggsAnalysisCombinedLimit_INCLUDE_DIRS
    "${HiggsAnalysisCombinedLimit_PREFIX}/include")
  set(HiggsAnalysisCombinedLimit_LIBRARIES
    "${HiggsAnalysisCombinedLimit_LIBRARY}")
endif()

mark_as_advanced(
  HiggsAnalysisCombinedLimit_INCLUDE_DIR
  HiggsAnalysisCombinedLimit_LIBRARY
  HiggsAnalysisCombinedLimit_PCM
  HiggsAnalysisCombinedLimit_ROOTMAP)
