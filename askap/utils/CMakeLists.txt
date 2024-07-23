add_sources_to_yandasoft(
	CleanUtils.cc
	LinmosUtils.cc
    EigenSolve.cc
	IlluminationUtils.cc
	ImplCalWeightSolver.cc
	SkyCatalogTabWriter.cc
    DelaySolverImpl.cc
    StatsAndMask.cc
	TilingUtils.cc
)

install (FILES
	CleanUtils.h
	LinmosUtils.h
	EigenSolve.h
	IlluminationUtils.h
	ImplCalWeightSolver.h
	SkyCatalogTabWriter.h
    DelaySolverImpl.h
    StatsAndMask.h
	TilingUtils.h
    DESTINATION include/askap/utils
)
