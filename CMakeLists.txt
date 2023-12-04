cmake_minimum_required (VERSION 2.6)
project (BRIQX)

# The version number.
set (BRIQX_VERSION_MAJOR 1)
set (BRIQX_VERSION_MINOR 0)

SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -L ${LIBRARY_OUTPUT_PATH}" )


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_library(briqx
dataio/src/datapaths.cpp
forcefield/src/AtomicClashEnergy.cpp
forcefield/src/BasePair6DEnergyTableCG.cpp
forcefield/src/BasePair6DEnergyTable.cpp
forcefield/src/ForceFieldPara.cpp
forcefield/src/HbondEnergy.cpp
forcefield/src/PO3Builder.cpp
forcefield/src/ProAtomicEnergyTable.cpp
forcefield/src/ProEnergyTable.cpp
forcefield/src/ProS1S2Energy.cpp
forcefield/src/ProteinRNAEnergyTable.cpp
forcefield/src/RiboseOxygenEnergyTable.cpp
forcefield/src/RnaAtomicEnergyTable.cpp
forcefield/src/RnaEnergyTable.cpp
forcefield/src/RnaEnergyTableSimple.cpp
forcefield/src/XPara.cpp
forcefield/test/CalRNAEnergy.cpp
forcefield/test/TestAtomicEnergy.cpp
forcefield/test/TestBaseAtomNameAndSep.cpp
forcefield/test/TestBasePairEnergy.cpp
forcefield/test/TestProPolarEnergy.cpp
forcefield/test/TestS1S2.cpp
geometry/src/Angles.cpp
geometry/src/CsMove.cpp
geometry/src/localframe.cpp
geometry/src/RMSD.cpp
geometry/src/TransMatrix.cpp
geometry/test/TestCoordSystem.cpp
geometry/test/TestTime.cpp
math/src/AAProbabilityArray.cpp
math/src/AAProbabilityMatrix.cpp
math/src/AAScoreArray.cpp
math/src/AAScoreMatrix.cpp
math/src/KDTree.cpp
math/src/Stat.cpp
math/src/SubstitutionMatrix.cpp
math/test/test.cpp
model/src/AssignProSSAndSasa.cpp
model/src/AssignRNASasa.cpp
model/src/AssignRNASS.cpp
model/src/AtomLib.cpp
model/src/BaseDistanceMatrix.cpp
model/src/BasePairComposition.cpp
model/src/BasePair.cpp
model/src/BasePairLib.cpp
model/src/BaseRotamer.cpp
model/src/BaseRotamerLib.cpp
model/src/DistanceMatrixHbond.cpp
model/src/DistanceMatrixRes.cpp
model/src/PhosphateRotamer.cpp
model/src/PhosphateRotamerLib.cpp
model/src/ResBBRotamer.cpp
model/src/ResBBRotamerLib.cpp
model/src/ResConformer.cpp
model/src/ResName.cpp
model/src/ResPeptideRotamer.cpp
model/src/ResPeptideRotamerLib.cpp
model/src/ResScRotamer.cpp
model/src/ResScRotamerLib.cpp
model/src/ResScRotamerLibMini.cpp
model/src/RiboseRotamer.cpp
model/src/RiboseRotamerLib.cpp
model/src/RNABaseLib.cpp
model/src/RNABaseName.cpp
model/src/RNARibose.cpp
model/src/RotamerLib.cpp
model/src/StructureModel.cpp
model/test/TestAtomLib.cpp
model/test/TestBackboneRotamer.cpp
model/test/TestBaseLib.cpp
model/test/TestBasePairLib.cpp
model/test/TestDM.cpp
model/test/TestRiboseRotamer.cpp
model/test/TestRotLib.cpp
model/test/TestSSAndSasa.cpp
motif/src/MotifAlign.cpp
motif/src/RGEdge.cpp
motif/src/RGNode.cpp
motif/src/RNAGraph.cpp
motif/src/SimpleAlign.cpp
motif/test/CalculateAllToAllScore.cpp
motif/test/rnaGraph.cpp
motif/test/TestAlignment.cpp
motif/test/TestAlignScore.cpp
motif/test/TestSimpleAlign.cpp
motif/test/writeMotifPDB.cpp
motif/test/writeSingleMotifPDB.cpp
para/src/DesignPara.cpp
para/src/ProParameter.cpp
predNA/src/BackboneModelingTemplate.cpp
predNA/src/BaseMoveLibrary.cpp
predNA/src/BRConnection.cpp
predNA/src/BRFoldingTree.cpp
predNA/src/BRNode.cpp
predNA/src/CsMoveTo6DKey.cpp
predNA/src/EdgeInformation.cpp
predNA/src/FragmentAssembling.cpp
predNA/src/FragmentLibrary.cpp
predNA/src/MCRun.cpp
predNA/src/MoveMutator.cpp
predNA/src/NuEnergyCalculator.cpp
predNA/src/NuGraph.cpp
predNA/src/NuMoveSet.cpp
predNA/src/OrientationIndex.cpp
predNA/src/ThreeBaseMoveLibrary.cpp
predNA/test/AddBaseToTarget.cpp
predNA/test/AssembleMultipleFragments.cpp
predNA/test/AssembleTwoFragments.cpp
predNA/test/BuildBackbone.cpp
predNA/test/BuildFoldingTree.cpp
predNA/test/BuildMutation.cpp
predNA/test/ConnectInfoToRMS.cpp
predNA/test/energyDebug.cpp
predNA/test/ExtendHashkey.cpp
predNA/test/ExtendRiboConnectHashKey.cpp
predNA/test/GenerateDecoysRandInit.cpp
predNA/test/GenerateF2FragmentLib.cpp
predNA/test/GenerateF3FragmentLib.cpp
predNA/test/GenerateInitModel.cpp
predNA/test/KeyToCsMove2.cpp
predNA/test/KeyToCsMove.cpp
predNA/test/KeyToHbondCsMove.cpp
predNA/test/optModel.cpp
predNA/test/ReplaceFragment.cpp
predNA/test/ReplaceMotif.cpp
predNA/test/RNARMSD.cpp
predNA/test/RNAScoring.cpp
predNA/test/RunPredict.cpp
predNA/test/RunRefinement.cpp
predNA/test/RunRefinementFast.cpp
predNA/test/runRnaXp.cpp
predNA/test/SimpleMC.cpp
predNA/test/TestAssignment.cpp
predNA/test/TestBackbone.cpp
predNA/test/TestConnectionEnergy.cpp
predNA/test/TestFoldTree.cpp
predNA/test/TestFragment.cpp
predNA/test/TestGraph.cpp
predNA/test/TestKey.cpp
predNA/test/TestMCSimPara.cpp
predNA/test/TestMerge.cpp
predNA/test/TestScoring.cpp
predNA/test/TrainPhoRep.cpp
protein/src/BuildMutationTemplate.cpp
protein/src/EnergyArray.cpp
protein/src/EnergyCalculator.cpp
protein/src/EnergyMatrix.cpp
protein/src/LoopModelingTemplate.cpp
protein/src/ResConnection.cpp
protein/src/ResFoldingTree.cpp
protein/src/ResMutator.cpp
protein/src/ResNode.cpp
protein/src/RotSequence.cpp
protein/src/SeqDesignTemplate.cpp
protein/src/SeqDesignTemplateFast.cpp
protein/src/SidechainModelingTemplate.cpp
protein/src/SingleSidechainModelingTemplate.cpp
protein/test/CalculateDetailEnergy.cpp
protein/test/Design1ubq.cpp
protein/test/GenerateSM.cpp
protein/test/KeyToHbondCsMove.cpp
protein/test/SinglePackToPDB.cpp
protein/test/TestAtomicEnergy.cpp
protein/test/TestDesign.cpp
protein/test/TestLoopModeling.cpp
protein/test/TestScModeling.cpp
protein/test/TestSinglechainFP.cpp
protein/test/TestSinglechainFV.cpp
protein/test/TestSingleMutation.cpp
protein/test/TestSingleSidechain.cpp
protein/test/TetsDesignRanks.cpp
tools/src/CmdArgs.cpp
tools/src/InputParser.cpp
tools/src/StringTool.cpp
)


add_subdirectory (model/test)
add_subdirectory (forcefield/test)
add_subdirectory (protein/test)
add_subdirectory (predNA/test)
add_subdirectory (motif/test)
