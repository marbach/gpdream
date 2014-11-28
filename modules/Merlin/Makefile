LFLAG = -lgsl -lgslcblas 
SRC = src/Framework.C src/MetaLearner.C src/MetaMove.C src/MotifManager.C src/BFGSWrapper.C src/common/Error.C src/common/EvidenceManager.C src/common/FactorManager.C src/common/LatticeStructure.C src/common/Potential.C src/common/SlimFactor.C src/common/VariableManager.C src/common/Evidence.C src/common/FactorGraph.C src/common/Matrix.C src/common/PotentialManager.C src/common/Variable.C src/HierarchicalCluster.C src/HierarchicalClusterNode.C src/common/Graph.C src/common/Vertex.C src/HyperGeomPval.C src/Distance.C src/Heap.C
LIBPATH = ../../../thirdparty/gsl_install/lib
INCLPATH1 = ../../..//thirdparty/gsl_install/include
INCLPATH2 = src/common

CC=g++
CFLAGS = -g

merlin: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o merlin

