cmake_minimum_required(VERSION 3.20)

project(CombineHarvester VERSION 4.0.0 LANGUAGES CXX)

option(CH_INSTALL_PYTHON
  "Install the existing CombineHarvester Python interface" ON)
set(CH_PYTHON_EXECUTABLE "" CACHE FILEPATH
  "Python interpreter whose purelib directory receives the interface")
set(CH_PYTHON_INSTALL_DIR "" CACHE STRING
  "Relative Python install directory; derived from the selected interpreter when empty")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  message(FATAL_ERROR
    "An explicit CMAKE_INSTALL_PREFIX is required. Use a disposable staging "
    "prefix; do not install into CONDA_PREFIX by default.")
endif()

set(CH_CONDA_PREFIX "$ENV{CONDA_PREFIX}" CACHE PATH
  "Approved dependency prefix (normally the combine-11 conda environment)")
if(NOT IS_DIRECTORY "${CH_CONDA_PREFIX}")
  message(FATAL_ERROR
    "CH_CONDA_PREFIX='${CH_CONDA_PREFIX}' is not a directory. Configure from "
    "combine-11 or pass -DCH_CONDA_PREFIX=/home/mohammad/.conda/envs/combine-11.")
endif()

# The conda compiler activation injects its own absolute library RPATH into
# CMake's global linker flags.  That is useful for an in-environment build but
# would make the installed CH artifacts permanently point at this particular
# conda checkout.  External dependencies are supplied by the selected
# environment at build/run time; retain -L and -rpath-link, and remove only
# this exact injected RPATH.
set(_ch_conda_rpath_flag "-Wl,-rpath,${CH_CONDA_PREFIX}/lib")
foreach(_ch_linker_flag_var IN ITEMS
    CMAKE_EXE_LINKER_FLAGS
    CMAKE_SHARED_LINKER_FLAGS
    CMAKE_MODULE_LINKER_FLAGS)
  if(DEFINED ${_ch_linker_flag_var})
    string(REPLACE "${_ch_conda_rpath_flag}" ""
      _ch_sanitized_linker_flags "${${_ch_linker_flag_var}}")
    set(${_ch_linker_flag_var} "${_ch_sanitized_linker_flags}")
    if(DEFINED CACHE{${_ch_linker_flag_var}})
      set(${_ch_linker_flag_var} "${_ch_sanitized_linker_flags}" CACHE STRING
        "${_ch_linker_flag_var}" FORCE)
    endif()
  endif()
endforeach()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

# The combine-11 GCC driver has a conda-forge default RPATH compiled into its
# specs file.  Remove that environment-specific path from each CH artifact
# after linking, otherwise it can take precedence over a CMSSW runtime's
# selected libraries.  patchelf is used only on the locally produced ELF
# files; no dependency is copied or rebuilt.
find_program(CH_PATCHELF_EXECUTABLE NAMES patchelf REQUIRED)
message(STATUS "ELF RPATH normalizer: ${CH_PATCHELF_EXECUTABLE}")
function(ch_normalize_rpath target relative_rpath)
  add_custom_command(TARGET "${target}" POST_BUILD
    COMMAND "${CH_PATCHELF_EXECUTABLE}" --set-rpath "${relative_rpath}"
            "$<TARGET_FILE:${target}>"
    COMMENT "Normalizing ${target} RPATH to ${relative_rpath}"
    VERBATIM)
endfunction()

include(GNUInstallDirs)
list(PREPEND CMAKE_PREFIX_PATH "${CH_CONDA_PREFIX}")
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(ROOT CONFIG REQUIRED COMPONENTS
  Core RIO Hist Graf Gpad Tree Matrix MathCore MathMore
  RooFitCore RooFit RooStats HistFactory)
find_package(Boost CONFIG REQUIRED COMPONENTS filesystem program_options regex)
find_package(Eigen3 CONFIG REQUIRED NO_MODULE)

# FindLibXml2 is used deliberately, but its search is pinned to the approved
# prefix and checked below so a host installation cannot be mixed in silently.
set(LibXml2_ROOT "${CH_CONDA_PREFIX}" CACHE PATH
  "LibXml2 prefix used by this build")
find_package(LibXml2 REQUIRED)
foreach(_xml_path IN ITEMS LIBXML2_INCLUDE_DIR LIBXML2_LIBRARY)
  if(NOT DEFINED ${_xml_path} OR NOT IS_ABSOLUTE "${${_xml_path}}")
    message(FATAL_ERROR "LibXml2 did not expose an absolute ${_xml_path}")
  endif()
  string(FIND "${${_xml_path}}" "${CH_CONDA_PREFIX}" _xml_prefix_offset)
  if(NOT _xml_prefix_offset EQUAL 0)
    message(FATAL_ERROR
      "LibXml2 ${_xml_path}='${${_xml_path}}' is outside CH_CONDA_PREFIX='${CH_CONDA_PREFIX}'")
  endif()
endforeach()

set(HiggsAnalysisCombinedLimit_ROOT "${CH_CONDA_PREFIX}" CACHE PATH
  "Installed CombinedLimit prefix; this project never builds CombinedLimit")
find_package(HiggsAnalysisCombinedLimit REQUIRED)
find_package(VDT REQUIRED)

if(CH_INSTALL_PYTHON)
  if(CH_PYTHON_EXECUTABLE)
    set(Python3_EXECUTABLE "${CH_PYTHON_EXECUTABLE}" CACHE FILEPATH
      "Python interpreter selected for CombineHarvester" FORCE)
  endif()
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  get_filename_component(_ch_python_real_executable
    "${Python3_EXECUTABLE}" REALPATH)
  string(FIND "${_ch_python_real_executable}" "${CH_CONDA_PREFIX}"
    _ch_python_prefix_offset)
  if(NOT _ch_python_prefix_offset EQUAL 0)
    message(FATAL_ERROR
      "Python3_EXECUTABLE='${Python3_EXECUTABLE}' is outside "
      "CH_CONDA_PREFIX='${CH_CONDA_PREFIX}'. Pass "
      "-DCH_PYTHON_EXECUTABLE=${CH_CONDA_PREFIX}/bin/python.")
  endif()

  if(NOT CH_PYTHON_INSTALL_DIR)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c
        "import sysconfig; print(sysconfig.get_path('purelib'))"
      RESULT_VARIABLE _ch_python_site_result
      OUTPUT_VARIABLE _ch_python_site_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_VARIABLE _ch_python_site_error)
    if(NOT _ch_python_site_result EQUAL 0 OR NOT IS_ABSOLUTE
       "${_ch_python_site_dir}")
      message(FATAL_ERROR
        "Could not determine Python purelib directory with "
        "${Python3_EXECUTABLE}: ${_ch_python_site_error}")
    endif()
    file(RELATIVE_PATH CH_PYTHON_INSTALL_DIR "${CH_CONDA_PREFIX}"
      "${_ch_python_site_dir}")
  endif()
  if(IS_ABSOLUTE "${CH_PYTHON_INSTALL_DIR}" OR
     "${CH_PYTHON_INSTALL_DIR}" MATCHES "(^|/)\.\.(\/|$)")
    message(FATAL_ERROR
      "CH_PYTHON_INSTALL_DIR must be a relative path inside the selected "
      "conda prefix, got '${CH_PYTHON_INSTALL_DIR}'.")
  endif()
  message(STATUS "Python executable: ${Python3_EXECUTABLE}")
  message(STATUS "Python install directory: ${CH_PYTHON_INSTALL_DIR}")

  set(_ch_generated_python_root
    "${CMAKE_CURRENT_BINARY_DIR}/generated-python/CombineHarvester")
  file(MAKE_DIRECTORY
    "${_ch_generated_python_root}/CombineTools"
    "${_ch_generated_python_root}/CombineTools/combine"
    "${_ch_generated_python_root}/CombinePdfs")
  foreach(_ch_python_init IN ITEMS
      "${_ch_generated_python_root}/__init__.py"
      "${_ch_generated_python_root}/CombineTools/__init__.py"
      "${_ch_generated_python_root}/CombineTools/combine/__init__.py"
      "${_ch_generated_python_root}/CombinePdfs/__init__.py")
    file(WRITE "${_ch_python_init}"
      "# Generated package marker for the installed CombineHarvester interface.\n")
  endforeach()
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CombineHarvesterPythonCh.py.in"
    "${_ch_generated_python_root}/CombineTools/ch.py" @ONLY)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

message(STATUS "CombineHarvester source: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Approved dependency prefix: ${CH_CONDA_PREFIX}")
message(STATUS "ROOT version: ${ROOT_VERSION}")
message(STATUS "ROOT prefix: ${ROOT_DIR}")
message(STATUS "Boost prefix: ${Boost_DIR}")
message(STATUS "Eigen3 prefix: ${Eigen3_DIR}")
message(STATUS "LibXml2 include: ${LIBXML2_INCLUDE_DIR}")
message(STATUS "LibXml2 library: ${LIBXML2_LIBRARY}")
message(STATUS "VDT include: ${VDT_INCLUDE_DIR}")
message(STATUS "VDT library: ${VDT_LIBRARY}")
message(STATUS "CombinedLimit include: ${HiggsAnalysisCombinedLimit_INCLUDE_DIR}")
message(STATUS "CombinedLimit library: ${HiggsAnalysisCombinedLimit_LIBRARY}")

include(CTest)
add_subdirectory(CombineTools)

install(TARGETS CombineHarvesterCombineTools ChronoSpectra
  EXPORT CombineHarvesterTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(EXPORT CombineHarvesterTargets
  FILE CombineHarvesterTargets.cmake
  NAMESPACE CombineHarvester::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/CombineHarvester")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/CombineTools/interface/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CombineHarvester/CombineTools/interface"
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp")

if(CH_INSTALL_PYTHON)
  # The standalone ROOT dictionary retains these headers in its payload. Keep
  # them in the install tree so cppyy can resolve free-function declarations
  # without reaching back into the source or build directory.
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CombineTools/CombineHarvesterDictionarySelection.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CombineHarvester/CombineTools")
  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CombineTools/src/classes.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/CombineTools/src/classes_def.xml"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CombineHarvester/CombineTools/src")

  set(_ch_python_install_root
    "${CH_PYTHON_INSTALL_DIR}/CombineHarvester")
  foreach(_ch_python_component IN ITEMS CombineTools CombinePdfs)
    install(DIRECTORY
      "${CMAKE_CURRENT_SOURCE_DIR}/${_ch_python_component}/python/"
      DESTINATION "${_ch_python_install_root}/${_ch_python_component}"
      FILES_MATCHING
        PATTERN "*.py"
        PATTERN "__pycache__" EXCLUDE)
  endforeach()
  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CombineTools/python/ch.py"
    DESTINATION "${_ch_python_install_root}/CombineTools"
    RENAME "_ch_original.py")
  install(DIRECTORY "${_ch_generated_python_root}/"
    DESTINATION "${CH_PYTHON_INSTALL_DIR}/CombineHarvester"
    FILES_MATCHING PATTERN "*.py")
endif()

include(CMakePackageConfigHelpers)
set(_ch_package_dir "${CMAKE_INSTALL_LIBDIR}/cmake/CombineHarvester")
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CombineHarvesterConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CombineHarvesterConfig.cmake"
  INSTALL_DESTINATION "${_ch_package_dir}")
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CombineHarvesterConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/CombineHarvesterConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/CombineHarvesterConfigVersion.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindHiggsAnalysisCombinedLimit.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindVDT.cmake"
  DESTINATION "${_ch_package_dir}")

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/lib/libCombineHarvesterCombineTools_rdict.pcm"
  "${CMAKE_CURRENT_BINARY_DIR}/lib/libCombineHarvesterCombineTools.rootmap"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")
