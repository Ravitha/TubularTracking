cmake_minimum_required(VERSION 3.19)

project(MLTubularTracking)

find_package(MeVisLab COMPONENTS ML MLBase Tools VesselGraph HINTS $ENV{MLAB_ROOT} REQUIRED)
find_package(Boost REQUIRED)
find_package(FME COMPONENTS MLGaussFiltersStable REQUIRED)
find_package(newmat REQUIRED)

add_library(MLTubularTracking "")
add_library(FME::MLTubularTracking ALIAS MLTubularTracking)

target_include_directories(MLTubularTracking PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# set_target_properties(MLTubularTracking PROPERTIES MEVIS_MAINTAINER hmeine)  # not implemented
# set_target_properties(MLTubularTracking PROPERTIES WARN_LEVEL LOW)  # not implemented
target_sources(MLTubularTracking
  PRIVATE
    MLTubularTrackingInit.h
    MLTubularTrackingSystem.h
    mlTubularTracking.h
    VesselTracker.h
    VesselGraph.h
    VesselData.h
    VesselProfile.h
    VesselTemplate.h
    VesselTemplateTracker.h
    ImPatch.h
    MLTubularTrackingInit.cpp
    mlTubularTracking.cpp
    VesselTracker.cpp
    VesselData.cpp
    VesselTemplateTracker.cpp
    VesselTemplate.cpp
    VesselGraph.cpp
)

target_link_libraries(MLTubularTracking
  PRIVATE
    Boost::headers
    FME::MLGaussFiltersStable
    MeVisLab::ML
    MeVisLab::MLBase
    MeVisLab::Tools
    MeVisLab::VesselGraph
    newmat::newmat
)

target_compile_definitions(MLTubularTracking
  PRIVATE
    MEVIS_TARGET=MLTubularTracking
    MLTUBULARTRACKING_EXPORTS
)

# relatedFiles:
#   ../../../Modules/ML/MLTubularTracking/MLTubularTracking.def
mlab_install(MLTubularTracking NS FME)
mlab_install_headers(MLTubularTracking)