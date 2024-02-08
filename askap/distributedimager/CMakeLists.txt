add_sources_to_yandasoft(
AdviseDI.cc
CalcCore.cc
ContinuumImager.cc
ContinuumMaster.cc
ContinuumWorker.cc
CubeComms.cc
MSSplitter.cc
WorkUnitContainer.cc
DataSourceManager.cc
WorkUnitCache.cc
)

install (FILES
AdviseDI.h
CalcCore.h
ContinuumImager.h
ContinuumMaster.h
ContinuumWorker.h
CubeBuilder.h
CubeBuilder.tcc
CubeComms.h
CubeManager.h
MSSplitter.h
WorkUnitContainer.h
DataSourceManager.h
WorkUnitCache.h
DESTINATION include/askap/distributedimager
)
