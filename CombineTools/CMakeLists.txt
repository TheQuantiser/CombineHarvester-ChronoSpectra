set(_ch_combinetools_sources
  "${CMAKE_CURRENT_SOURCE_DIR}/src/AZhSystematics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/AutoRebin.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/BinByBin.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CardWriter.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CombineHarvester.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CombineHarvester_Creation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CombineHarvester_Datacards.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CombineHarvester_Evaluate.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CombineHarvester_Filters.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/CopyTools.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HhhSystematics.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HistMapping.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HttSystematics_MSSMLegacy.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HttSystematics_MSSMUpdate.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/HttSystematics_SMLegacy.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/JsonTools.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Logging.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Object.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Observation.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Parameter.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/ParseCombineWorkspace.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Process.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/SOverBTools.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Systematic.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/TFileIO.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Utilities.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationTools.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/json.cc")

set(_ch_generated_include "${CMAKE_CURRENT_BINARY_DIR}/generated-include")
file(MAKE_DIRECTORY "${_ch_generated_include}/CombineHarvester")
if(NOT EXISTS "${_ch_generated_include}/CombineHarvester/CombineTools")
  file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}"
    "${_ch_generated_include}/CombineHarvester/CombineTools" SYMBOLIC)
endif()

set(_ch_root_libraries
  ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree ROOT::Matrix
  ROOT::MathCore ROOT::MathMore ROOT::RooFitCore ROOT::RooFit ROOT::RooStats
  ROOT::HistFactory)

add_library(CombineHarvesterCombineTools SHARED ${_ch_combinetools_sources})
add_library(CombineHarvester::CombineTools ALIAS CombineHarvesterCombineTools)
set_target_properties(CombineHarvesterCombineTools PROPERTIES
  EXPORT_NAME CombineTools
  OUTPUT_NAME CombineHarvesterCombineTools
  INSTALL_RPATH "$ORIGIN"
  INSTALL_RPATH_USE_LINK_PATH FALSE)
ch_normalize_rpath(CombineHarvesterCombineTools "$ORIGIN")
target_include_directories(CombineHarvesterCombineTools
  PUBLIC
    "$<BUILD_INTERFACE:${_ch_generated_include}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(CombineHarvesterCombineTools
  PUBLIC
    HiggsAnalysisCombinedLimit::HiggsAnalysisCombinedLimit
    VDT::VDT
    LibXml2::LibXml2
    Boost::filesystem
    Boost::regex
    ${_ch_root_libraries})

# ROOT 6.36 recursively emits a collection proxy for the private
# AutoMCStatsSettings map even though ch::CombineHarvester is transient. The
# wrapper history shows that moving the nested struct to public access makes
# rootcling proceed, but that would alter the CMSSW-facing source API. ROOT's
# supported ROOT::Meta::Selection metadata instead allows this one member's
# type to be excluded from automatic dictionary selection while retaining the
# outer class dictionary. Generate that metadata only in the build tree; the
# current classes.h and classes_def.xml remain the exact dictionary inputs.
set(_ch_dictionary_selection
  "${CMAKE_CURRENT_BINARY_DIR}/CombineHarvesterDictionarySelection.h")
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/CombineHarvesterDictionarySelection.h.in"
  "${_ch_dictionary_selection}" COPYONLY)
message(STATUS
  "ROOT dictionary uses build-tree-only kNoAutoSelected metadata for "
  "ch::CombineHarvester::auto_stats_settings_")

include("${ROOT_DIR}/RootMacros.cmake")
ROOT_GENERATE_DICTIONARY(G__CombineHarvesterCombineTools
  "${_ch_dictionary_selection}"
  CombineHarvester/CombineTools/src/classes.h
  LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml"
  MODULE CombineHarvesterCombineTools
  NOINSTALL
  OPTIONS --noIncludePaths)

set(_ch_dictionary_pcm
  "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libCombineHarvesterCombineTools_rdict.pcm")
set(_ch_dictionary_rootmap
  "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libCombineHarvesterCombineTools.rootmap")
set_source_files_properties("${_ch_dictionary_pcm}" "${_ch_dictionary_rootmap}"
  PROPERTIES GENERATED TRUE)

add_executable(ChronoSpectra "${CMAKE_CURRENT_SOURCE_DIR}/bin/ChronoSpectra.cpp")
set_target_properties(ChronoSpectra PROPERTIES
  INSTALL_RPATH "$ORIGIN/../lib"
  INSTALL_RPATH_USE_LINK_PATH FALSE)
ch_normalize_rpath(ChronoSpectra "$ORIGIN/../lib")
target_include_directories(ChronoSpectra PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/bin")
target_link_libraries(ChronoSpectra PRIVATE
  CombineHarvester::CombineTools
  Boost::program_options
  Boost::filesystem
  Boost::regex
  ${_ch_root_libraries})

add_executable(ChronoSpectraInternalTest
  "${CMAKE_CURRENT_SOURCE_DIR}/test/ChronoSpectraInternalTest.cpp")
target_include_directories(ChronoSpectraInternalTest PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/bin")
target_link_libraries(ChronoSpectraInternalTest PRIVATE
  CombineHarvester::CombineTools)
ch_normalize_rpath(ChronoSpectraInternalTest "$ORIGIN/../lib")

if(BUILD_TESTING)
  add_test(NAME ChronoSpectraInternalTest COMMAND ChronoSpectraInternalTest)
  set_tests_properties(ChronoSpectraInternalTest PROPERTIES
    ENVIRONMENT
      "LD_LIBRARY_PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:${CH_CONDA_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}")
endif()
