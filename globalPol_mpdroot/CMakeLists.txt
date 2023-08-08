set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdPid
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder/mlem
  ${CMAKE_SOURCE_DIR}/detectors/zdc
  ${CMAKE_SOURCE_DIR}/detectors/tof
  ${CMAKE_SOURCE_DIR}/physics
  ${CMAKE_SOURCE_DIR}/physics/globalPolarization
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
  MpdAnalysisTask2.cxx
  MpdLambdaPol.cxx
  MpdGlobalPolarizationMC.cxx
  MpdGlobalPolarizationRECO.cxx
)

Set(HEADERS)
Set(LINKDEF MpdGlobalPolarizationLinkDef.h)
Set(LIBRARY_NAME MpdGlobalPolarization)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
