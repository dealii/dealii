## ---------------------------------------------------------------------
##
## Copyright (C) 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# avoid some annoying warnings...
#

STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wall")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-pedantic")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wunused-parameter")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wfloat-conversion")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wsign-compare")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wformat")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Woverloaded-virtual")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wformat-invalid-specifier")


#
# Add header and source files
#

SET(hdr_mesquite
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_version.h
  src/Control/Instruction.hpp
  src/Control/InstructionQueue.hpp
  src/Control/IQInterface.hpp
  src/Control/Settings.hpp
  src/Control/SlaveBoundaryVertices.hpp
  src/Control/TerminationCriterion.hpp
  src/Control/VertexSlaver.hpp
  src/include/CurveDomain.hpp
  src/include/MeshInterface.hpp
  src/include/Mesquite.hpp
  src/include/ParallelHelperInterface.hpp
  src/include/ParallelMeshInterface.hpp
  src/MappingFunction/Lagrange/HexLagrangeShape.hpp
  src/MappingFunction/Lagrange/QuadLagrangeShape.hpp
  src/MappingFunction/Lagrange/TetLagrangeShape.hpp
  src/MappingFunction/Lagrange/TriLagrangeShape.hpp
  src/MappingFunction/Linear/LinearHexahedron.hpp
  src/MappingFunction/Linear/LinearPrism.hpp
  src/MappingFunction/Linear/LinearPyramid.hpp
  src/MappingFunction/Linear/LinearQuadrilateral.hpp
  src/MappingFunction/Linear/LinearTetrahedron.hpp
  src/MappingFunction/Linear/LinearTriangle.hpp
  src/MappingFunction/MappingFunction.hpp
  src/MappingFunction/NodeSet.hpp
  src/MappingFunction/Sample.hpp
  src/Mesh/ArrayMesh.hpp
  src/Mesh/ElementPatches.hpp
  src/Mesh/ExtraData.hpp
  src/Mesh/ExtraDataUser.hpp
  src/Mesh/GlobalPatch.hpp
  src/Mesh/IdealElements.hpp
  src/Mesh/MeshDecorator.hpp
  src/Mesh/MeshImpl.hpp
  src/Mesh/MeshImplData.hpp
  src/Mesh/MeshImplTags.hpp
  src/Mesh/MeshUtil.hpp
  src/Mesh/MeshWriter.hpp
  src/Mesh/MsqFreeVertexIndexIterator.hpp
  src/Mesh/MsqMeshEntity.hpp
  src/Mesh/MsqVertex.hpp
  # src/Mesh/ParallelHelper.hpp
  # src/Mesh/ParallelHelperNoMPIStub.hpp
  src/Mesh/ParallelMeshImpl.hpp
  src/Mesh/PatchData.hpp
  src/Mesh/PatchIterator.hpp
  src/Mesh/PatchSet.hpp
  src/Mesh/TagVertexMesh.hpp
  src/Mesh/TopologyInfo.hpp
  src/Mesh/VertexPatches.hpp
  src/Misc/Bits.hpp
  src/Misc/BoundedCylinderDomain.hpp
  src/Misc/ConicDomain.hpp
  src/Misc/CylinderDomain.hpp
  src/Misc/DomainClassifier.hpp
  src/Misc/DomainUtil.hpp
  src/Misc/EdgeIterator.hpp
  src/Misc/Exponent.hpp
  src/Misc/FileTokenizer.hpp
  src/Misc/Matrix3D.hpp
  src/Misc/MeshDomain1D.hpp
  src/Misc/MeshTransform.hpp
  src/Misc/MsqDebug.hpp
  src/Misc/MsqError.hpp
  src/Misc/MsqFPE.hpp
  src/Misc/MsqGeomPrim.hpp
  src/Misc/MsqHessian.hpp
  # src/Misc/MsqIBase.hpp
  src/Misc/MsqInterrupt.hpp
  src/Misc/MsqMatrix.hpp
  src/Misc/MsqTimer.hpp
  src/Misc/PlanarDomain.hpp
  src/Misc/SimpleStats.hpp
  src/Misc/SphericalDomain.hpp
  src/Misc/SymMatrix3D.hpp
  src/Misc/Vector3D.hpp
  src/Misc/VtkTypeInfo.hpp
  src/Misc/XYPlanarDomain.hpp
  src/Misc/XYRectangle.hpp
  src/ObjectiveFunction/CompositeOFAdd.hpp
  src/ObjectiveFunction/CompositeOFMultiply.hpp
  src/ObjectiveFunction/CompositeOFScalarAdd.hpp
  src/ObjectiveFunction/CompositeOFScalarMultiply.hpp
  src/ObjectiveFunction/LInfTemplate.hpp
  src/ObjectiveFunction/LPtoPTemplate.hpp
  src/ObjectiveFunction/MaxTemplate.hpp
  src/ObjectiveFunction/ObjectiveFunction.hpp
  src/ObjectiveFunction/ObjectiveFunctionTemplate.hpp
  src/ObjectiveFunction/OFEvaluator.hpp
  src/ObjectiveFunction/PatchPowerMeanP.hpp
  src/ObjectiveFunction/PMeanPTemplate.hpp
  src/ObjectiveFunction/StdDevTemplate.hpp
  src/ObjectiveFunction/VarianceTemplate.hpp
  src/QualityAssessor/QualityAssessor.hpp
  src/QualityImprover/NullImprover.hpp
  src/QualityImprover/OptSolvers/ConjugateGradient.hpp
  src/QualityImprover/OptSolvers/FeasibleNewton.hpp
  src/QualityImprover/OptSolvers/NonGradient.hpp
  src/QualityImprover/OptSolvers/NonSmoothDescent.hpp
  src/QualityImprover/OptSolvers/QuasiNewton.hpp
  src/QualityImprover/OptSolvers/SteepestDescent.hpp
  src/QualityImprover/OptSolvers/TrustRegion.hpp
  src/QualityImprover/PatchSetUser.hpp
  src/QualityImprover/QualityImprover.hpp
  src/QualityImprover/Relaxation/LaplacianSmoother.hpp
  src/QualityImprover/Relaxation/Randomize.hpp
  src/QualityImprover/Relaxation/RelaxationSmoother.hpp
  src/QualityImprover/Relaxation/SmartLaplacianSmoother.hpp
  src/QualityImprover/VertexMover.hpp
  src/QualityMetric/AddQualityMetric.hpp
  src/QualityMetric/AveragingQM.hpp
  src/QualityMetric/Debug/CompareQM.hpp
  src/QualityMetric/Debug/NumericalQM.hpp
  src/QualityMetric/EdgeQM.hpp
  src/QualityMetric/ElementAvgQM.hpp
  src/QualityMetric/ElementMaxQM.hpp
  src/QualityMetric/ElementPMeanP.hpp
  src/QualityMetric/ElementQM.hpp
  src/QualityMetric/ElemSampleQM.hpp
  src/QualityMetric/MultiplyQualityMetric.hpp
  src/QualityMetric/PMeanPMetric.hpp
  src/QualityMetric/PowerQualityMetric.hpp
  src/QualityMetric/QualityMetric.hpp
  src/QualityMetric/ScalarAddQualityMetric.hpp
  src/QualityMetric/ScalarMultiplyQualityMetric.hpp
  src/QualityMetric/Shape/AspectRatioGammaQualityMetric.hpp
  src/QualityMetric/Shape/ConditionNumberFunctions.hpp
  src/QualityMetric/Shape/ConditionNumberQualityMetric.hpp
  src/QualityMetric/Shape/IdealWeightInverseMeanRatio.hpp
  src/QualityMetric/Shape/IdealWeightMeanRatio.hpp
  src/QualityMetric/Shape/MeanRatioFunctions.hpp
  src/QualityMetric/Shape/VertexConditionNumberQualityMetric.hpp
  src/QualityMetric/Smoothness/EdgeLengthQualityMetric.hpp
  src/QualityMetric/Smoothness/EdgeLengthRangeQualityMetric.hpp
  src/QualityMetric/TMP/AffineMapMetric.hpp
  src/QualityMetric/TMP/AWQualityMetric.hpp
  src/QualityMetric/TMP/TargetMetricUtil.hpp
  src/QualityMetric/TMP/TMPQualityMetric.hpp
  src/QualityMetric/TMP/TQualityMetric.hpp
  src/QualityMetric/Untangle/UntangleBetaQualityMetric.hpp
  src/QualityMetric/VertexMaxQM.hpp
  src/QualityMetric/VertexPMeanP.hpp
  src/QualityMetric/VertexQM.hpp
  src/QualityMetric/Volume/EdgeLengthMetric.hpp
  src/QualityMetric/Volume/LocalSizeQualityMetric.hpp
  src/QualityMetric/Volume/SizeMetric.hpp
  src/TargetCalculator/CachingTargetCalculator.hpp
  src/TargetCalculator/IdealShapeTarget.hpp
  src/TargetCalculator/InverseMetricWeight.hpp
  src/TargetCalculator/JacobianCalculator.hpp
  src/TargetCalculator/LambdaConstant.hpp
  src/TargetCalculator/LambdaTarget.hpp
  src/TargetCalculator/LVQDTargetCalculator.hpp
  src/TargetCalculator/MetricWeight.hpp
  src/TargetCalculator/ReferenceMesh.hpp
  src/TargetCalculator/RefMeshTargetCalculator.hpp
  src/TargetCalculator/RefSizeTargetCalculator.hpp
  src/TargetCalculator/RemainingWeight.hpp
  src/TargetCalculator/TargetCalculator.hpp
  src/TargetCalculator/TargetReader.hpp
  src/TargetCalculator/TargetWriter.hpp
  src/TargetCalculator/TetDihedralWeight.hpp
  src/TargetCalculator/WeightCalculator.hpp
  src/TargetCalculator/WeightReader.hpp
  src/TargetMetric/AWMetric.hpp
  src/TargetMetric/AWMetricBarrier.hpp
  src/TargetMetric/AWMetricNonBarrier.hpp
  src/TargetMetric/Misc/InvTransBarrier.hpp
  src/TargetMetric/Misc/TMixed.hpp
  src/TargetMetric/Misc/TOffset.hpp
  src/TargetMetric/Misc/TPower2.hpp
  src/TargetMetric/Misc/TScale.hpp
  src/TargetMetric/Misc/TSquared.hpp
  src/TargetMetric/Misc/TSum.hpp
  src/TargetMetric/Shape/AWShape2DB1.hpp
  src/TargetMetric/Shape/AWShape2DNB1.hpp
  src/TargetMetric/Shape/AWShape2DNB2.hpp
  src/TargetMetric/Shape/TInverseMeanRatio.hpp
  src/TargetMetric/Shape/TShape2DNB2.hpp
  src/TargetMetric/Shape/TShape3DB2.hpp
  src/TargetMetric/Shape/TShapeB1.hpp
  src/TargetMetric/Shape/TShapeNB1.hpp
  src/TargetMetric/ShapeOrient/AWShapeOrientNB1.hpp
  src/TargetMetric/ShapeOrient/TShapeOrientB1.hpp
  src/TargetMetric/ShapeOrient/TShapeOrientB2.hpp
  src/TargetMetric/ShapeOrient/TShapeOrientNB1.hpp
  src/TargetMetric/ShapeOrient/TShapeOrientNB2.hpp
  src/TargetMetric/ShapeSize/AWShapeSizeB1.hpp
  src/TargetMetric/ShapeSize/TShapeSize2DB2.hpp
  src/TargetMetric/ShapeSize/TShapeSize2DNB1.hpp
  src/TargetMetric/ShapeSize/TShapeSize2DNB2.hpp
  src/TargetMetric/ShapeSize/TShapeSize3DB2.hpp
  src/TargetMetric/ShapeSize/TShapeSize3DB4.hpp
  src/TargetMetric/ShapeSize/TShapeSize3DNB1.hpp
  src/TargetMetric/ShapeSize/TShapeSizeB1.hpp
  src/TargetMetric/ShapeSize/TShapeSizeB3.hpp
  src/TargetMetric/ShapeSize/TShapeSizeNB3.hpp
  src/TargetMetric/ShapeSizeOrient/AWShapeSizeOrientNB1.hpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientB1.hpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientB2.hpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientNB1.hpp
  src/TargetMetric/Size/AWSizeB1.hpp
  src/TargetMetric/Size/AWSizeNB1.hpp
  src/TargetMetric/Size/TSizeB1.hpp
  src/TargetMetric/Size/TSizeNB1.hpp
  src/TargetMetric/Size/TTau.hpp
  src/TargetMetric/TMetric.hpp
  src/TargetMetric/TMetricBarrier.hpp
  src/TargetMetric/TMetricNonBarrier.hpp
  src/TargetMetric/TMPCommon.hpp
  src/TargetMetric/TMPDerivs.hpp
  src/TargetMetric/Untangle/AWUntangleBeta.hpp
  src/TargetMetric/Untangle/TUntangle1.hpp
  src/TargetMetric/Untangle/TUntangleBeta.hpp
  src/TargetMetric/Untangle/TUntangleMu.hpp
  src/Wrappers/DeformingDomainWrapper.hpp
  src/Wrappers/LaplaceWrapper.hpp
  src/Wrappers/PaverMinEdgeLengthWrapper.hpp
  src/Wrappers/ShapeImprovementWrapper.hpp
  src/Wrappers/ShapeImprover.hpp
  src/Wrappers/SizeAdaptShapeWrapper.hpp
  src/Wrappers/UntangleWrapper.hpp
  src/Wrappers/ViscousCFDTetShapeWrapper.hpp
  src/Wrappers/Wrapper.hpp
)

SET(src_mesquite
  src/Control/Instruction.cpp
  src/Control/InstructionQueue.cpp
  src/Control/IQInterface.cpp
  src/Control/Settings.cpp
  src/Control/SlaveBoundaryVertices.cpp
  src/Control/TerminationCriterion.cpp
  src/MappingFunction/Lagrange/HexLagrangeShape.cpp
  src/MappingFunction/Lagrange/QuadLagrangeShape.cpp
  src/MappingFunction/Lagrange/TetLagrangeShape.cpp
  src/MappingFunction/Lagrange/TriLagrangeShape.cpp
  src/MappingFunction/Linear/LinearHexahedron.cpp
  src/MappingFunction/Linear/LinearPrism.cpp
  src/MappingFunction/Linear/LinearPyramid.cpp
  src/MappingFunction/Linear/LinearQuadrilateral.cpp
  src/MappingFunction/Linear/LinearTetrahedron.cpp
  src/MappingFunction/Linear/LinearTriangle.cpp
  src/MappingFunction/MappingFunction.cpp
  src/MappingFunction/NodeSet.cpp
  src/MappingFunction/Sample.cpp
  src/Mesh/ArrayMesh.cpp
  src/Mesh/ElementPatches.cpp
  src/Mesh/ExtraData.cpp
  src/Mesh/GlobalPatch.cpp
  src/Mesh/IdealElements.cpp
  src/Mesh/MeshDecorator.cpp
  src/Mesh/MeshImpl.cpp
  src/Mesh/MeshImplData.cpp
  src/Mesh/MeshImplTags.cpp
  src/Mesh/MeshUtil.cpp
  src/Mesh/MeshWriter.cpp
  src/Mesh/MsqMeshEntity.cpp
  src/Mesh/MsqVertex.cpp
  # src/Mesh/ParallelHelper.cpp
  # src/Mesh/ParallelHelperNoMPIStub.cpp
  src/Mesh/ParallelMeshImpl.cpp
  src/Mesh/PatchData.cpp
  src/Mesh/PatchIterator.cpp
  src/Mesh/PatchSet.cpp
  src/Mesh/TagVertexMesh.cpp
  src/Mesh/TopologyInfo.cpp
  src/Mesh/VertexPatches.cpp
  src/Misc/BoundedCylinderDomain.cpp
  src/Misc/ConicDomain.cpp
  src/Misc/CylinderDomain.cpp
  src/Misc/DomainClassifier.cpp
  src/Misc/DomainUtil.cpp
  src/Misc/EdgeIterator.cpp
  src/Misc/Exponent.cpp
  src/Misc/FileTokenizer.cpp
  src/Misc/MeshDomain1D.cpp
  src/Misc/MeshTransform.cpp
  src/Misc/MesquiteVersion.cpp
  src/Misc/MsqDebug.cpp
  src/Misc/MsqError.cpp
  src/Misc/MsqFPE.cpp
  src/Misc/MsqGeomPrim.cpp
  src/Misc/MsqHessian.cpp
  # src/Misc/MsqIBase.cpp
  src/Misc/MsqInterrupt.cpp
  src/Misc/MsqTimer.cpp
  src/Misc/PlanarDomain.cpp
  src/Misc/SimpleStats.cpp
  src/Misc/SphericalDomain.cpp
  src/Misc/Vector3D.cpp
  src/Misc/VtkTypeInfo.cpp
  src/Misc/XYPlanarDomain.cpp
  src/Misc/XYRectangle.cpp
  src/ObjectiveFunction/CompositeOFAdd.cpp
  src/ObjectiveFunction/CompositeOFMultiply.cpp
  src/ObjectiveFunction/CompositeOFScalarAdd.cpp
  src/ObjectiveFunction/CompositeOFScalarMultiply.cpp
  src/ObjectiveFunction/LInfTemplate.cpp
  src/ObjectiveFunction/LPtoPTemplate.cpp
  src/ObjectiveFunction/MaxTemplate.cpp
  src/ObjectiveFunction/ObjectiveFunction.cpp
  src/ObjectiveFunction/ObjectiveFunctionTemplate.cpp
  src/ObjectiveFunction/OFEvaluator.cpp
  src/ObjectiveFunction/PatchPowerMeanP.cpp
  src/ObjectiveFunction/PMeanPTemplate.cpp
  src/ObjectiveFunction/StdDevTemplate.cpp
  src/ObjectiveFunction/VarianceTemplate.cpp
  src/QualityAssessor/QualityAssessor.cpp
  src/QualityImprover/OptSolvers/ConjugateGradient.cpp
  src/QualityImprover/OptSolvers/FeasibleNewton.cpp
  src/QualityImprover/OptSolvers/NonGradient.cpp
  src/QualityImprover/OptSolvers/NonSmoothDescent.cpp
  src/QualityImprover/OptSolvers/QuasiNewton.cpp
  src/QualityImprover/OptSolvers/SteepestDescent.cpp
  src/QualityImprover/OptSolvers/TrustRegion.cpp
  src/QualityImprover/PatchSetUser.cpp
  src/QualityImprover/QualityImprover.cpp
  src/QualityImprover/Relaxation/LaplacianSmoother.cpp
  src/QualityImprover/Relaxation/Randomize.cpp
  src/QualityImprover/Relaxation/RelaxationSmoother.cpp
  src/QualityImprover/Relaxation/SmartLaplacianSmoother.cpp
  src/QualityImprover/VertexMover.cpp
  src/QualityMetric/AddQualityMetric.cpp
  src/QualityMetric/AveragingQM.cpp
  src/QualityMetric/Debug/CompareQM.cpp
  src/QualityMetric/Debug/NumericalQM.cpp
  src/QualityMetric/EdgeQM.cpp
  src/QualityMetric/ElementAvgQM.cpp
  src/QualityMetric/ElementMaxQM.cpp
  src/QualityMetric/ElementPMeanP.cpp
  src/QualityMetric/ElementQM.cpp
  src/QualityMetric/ElemSampleQM.cpp
  src/QualityMetric/MultiplyQualityMetric.cpp
  src/QualityMetric/PMeanPMetric.cpp
  src/QualityMetric/PowerQualityMetric.cpp
  src/QualityMetric/QualityMetric.cpp
  src/QualityMetric/ScalarAddQualityMetric.cpp
  src/QualityMetric/ScalarMultiplyQualityMetric.cpp
  src/QualityMetric/Shape/AspectRatioGammaQualityMetric.cpp
  src/QualityMetric/Shape/ConditionNumberQualityMetric.cpp
  src/QualityMetric/Shape/IdealWeightInverseMeanRatio.cpp
  src/QualityMetric/Shape/IdealWeightMeanRatio.cpp
  src/QualityMetric/Shape/VertexConditionNumberQualityMetric.cpp
  src/QualityMetric/Smoothness/EdgeLengthQualityMetric.cpp
  src/QualityMetric/Smoothness/EdgeLengthRangeQualityMetric.cpp
  src/QualityMetric/TMP/AffineMapMetric.cpp
  src/QualityMetric/TMP/AWQualityMetric.cpp
  src/QualityMetric/TMP/TargetMetricUtil.cpp
  src/QualityMetric/TMP/TMPQualityMetric.cpp
  src/QualityMetric/TMP/TQualityMetric.cpp
  src/QualityMetric/Untangle/UntangleBetaQualityMetric.cpp
  src/QualityMetric/VertexMaxQM.cpp
  src/QualityMetric/VertexPMeanP.cpp
  src/QualityMetric/VertexQM.cpp
  src/QualityMetric/Volume/EdgeLengthMetric.cpp
  src/QualityMetric/Volume/LocalSizeQualityMetric.cpp
  src/QualityMetric/Volume/SizeMetric.cpp
  src/TargetCalculator/CachingTargetCalculator.cpp
  src/TargetCalculator/IdealShapeTarget.cpp
  src/TargetCalculator/InverseMetricWeight.cpp
  src/TargetCalculator/JacobianCalculator.cpp
  src/TargetCalculator/LambdaConstant.cpp
  src/TargetCalculator/LambdaTarget.cpp
  src/TargetCalculator/LVQDTargetCalculator.cpp
  src/TargetCalculator/MetricWeight.cpp
  src/TargetCalculator/ReferenceMesh.cpp
  src/TargetCalculator/RefMeshTargetCalculator.cpp
  src/TargetCalculator/RefSizeTargetCalculator.cpp
  src/TargetCalculator/RemainingWeight.cpp
  src/TargetCalculator/TargetCalculator.cpp
  src/TargetCalculator/TargetReader.cpp
  src/TargetCalculator/TargetWriter.cpp
  src/TargetCalculator/TetDihedralWeight.cpp
  src/TargetCalculator/WeightCalculator.cpp
  src/TargetCalculator/WeightReader.cpp
  src/TargetMetric/AWMetric.cpp
  src/TargetMetric/AWMetricBarrier.cpp
  src/TargetMetric/AWMetricNonBarrier.cpp
  src/TargetMetric/Misc/InvTransBarrier.cpp
  src/TargetMetric/Misc/TMixed.cpp
  src/TargetMetric/Misc/TOffset.cpp
  src/TargetMetric/Misc/TPower2.cpp
  src/TargetMetric/Misc/TScale.cpp
  src/TargetMetric/Misc/TSquared.cpp
  src/TargetMetric/Misc/TSum.cpp
  src/TargetMetric/Shape/AWShape2DB1.cpp
  src/TargetMetric/Shape/AWShape2DNB1.cpp
  src/TargetMetric/Shape/AWShape2DNB2.cpp
  src/TargetMetric/Shape/TInverseMeanRatio.cpp
  src/TargetMetric/Shape/TShape2DNB2.cpp
  src/TargetMetric/Shape/TShape3DB2.cpp
  src/TargetMetric/Shape/TShapeB1.cpp
  src/TargetMetric/Shape/TShapeNB1.cpp
  src/TargetMetric/ShapeOrient/AWShapeOrientNB1.cpp
  src/TargetMetric/ShapeOrient/TShapeOrientB1.cpp
  src/TargetMetric/ShapeOrient/TShapeOrientB2.cpp
  src/TargetMetric/ShapeOrient/TShapeOrientNB1.cpp
  src/TargetMetric/ShapeOrient/TShapeOrientNB2.cpp
  src/TargetMetric/ShapeSize/AWShapeSizeB1.cpp
  src/TargetMetric/ShapeSize/TShapeSize2DB2.cpp
  src/TargetMetric/ShapeSize/TShapeSize2DNB1.cpp
  src/TargetMetric/ShapeSize/TShapeSize2DNB2.cpp
  src/TargetMetric/ShapeSize/TShapeSize3DB2.cpp
  src/TargetMetric/ShapeSize/TShapeSize3DB4.cpp
  src/TargetMetric/ShapeSize/TShapeSize3DNB1.cpp
  src/TargetMetric/ShapeSize/TShapeSizeB1.cpp
  src/TargetMetric/ShapeSize/TShapeSizeB3.cpp
  src/TargetMetric/ShapeSize/TShapeSizeNB3.cpp
  src/TargetMetric/ShapeSizeOrient/AWShapeSizeOrientNB1.cpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientB1.cpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientB2.cpp
  src/TargetMetric/ShapeSizeOrient/TShapeSizeOrientNB1.cpp
  src/TargetMetric/Size/AWSizeB1.cpp
  src/TargetMetric/Size/AWSizeNB1.cpp
  src/TargetMetric/Size/TSizeB1.cpp
  src/TargetMetric/Size/TSizeNB1.cpp
  src/TargetMetric/Size/TTau.cpp
  src/TargetMetric/TMetric.cpp
  src/TargetMetric/TMetricBarrier.cpp
  src/TargetMetric/TMetricNonBarrier.cpp
  src/TargetMetric/Untangle/AWUntangleBeta.cpp
  src/TargetMetric/Untangle/TUntangle1.cpp
  src/TargetMetric/Untangle/TUntangleBeta.cpp
  src/TargetMetric/Untangle/TUntangleMu.cpp
  src/Wrappers/DeformingDomainWrapper.cpp
  src/Wrappers/LaplaceWrapper.cpp
  src/Wrappers/PaverMinEdgeLengthWrapper.cpp
  src/Wrappers/ShapeImprovementWrapper.cpp
  src/Wrappers/ShapeImprover.cpp
  src/Wrappers/SizeAdaptShapeWrapper.cpp
  src/Wrappers/UntangleWrapper.cpp
  src/Wrappers/ViscousCFDTetShapeWrapper.cpp
  src/Wrappers/Wrapper.cpp
)

IF (DEAL_II_WITH_MPI)
  MESSAGE(STATUS "Mesquite will be built with MPI support.")
  SET(hdr_mesquite ${hdr_mesquite} src/Mesh/ParallelHelper.hpp)
  SET(src_mesquite ${src_mesquite} src/Mesh/ParallelHelper.cpp)
ELSE()
  MESSAGE(STATUS "Mesquite will not be built with MPI support.")
  SET(hdr_mesquite ${hdr_mesquite} src/Mesh/ParallelHelperNoMPIStub.hpp)
  SET(src_mesquite ${src_mesquite} src/Mesh/ParallelHelperNoMPIStub.cpp)
ENDIF()

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Control
  ${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction
  ${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction/Linear
  ${CMAKE_CURRENT_SOURCE_DIR}/src/MappingFunction/Lagrange
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Misc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ObjectiveFunction
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityAssessor
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover/OptSolvers
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityImprover/Relaxation
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Debug
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Shape
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Smoothness
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Volume
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/Untangle
  ${CMAKE_CURRENT_SOURCE_DIR}/src/QualityMetric/TMP
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetCalculator
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Misc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Shape
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeOrient
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeSize
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/ShapeSizeOrient
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Size
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TargetMetric/Untangle
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Wrappers
)


#
# Create the Mesquite_all_headers.hpp file.
#   - Sourced from <mesquite>/src/CMakeLists.txt
#

SET( all_file ${CMAKE_CURRENT_BINARY_DIR}/include/Mesquite_all_headers.hpp )
# Cannot do this correctly because list if all headers is too long
# of a command-line argument for some platforms (e.g. Windows.)
# Instead, use less optimial solution (doesn't automatically update.)
#
#ADD_CUSTOM_COMMAND( 
#  OUTPUT ${all_file}
#  COMMAND ${CMAKE_COMMAND}
#  ARGS -D output_file:STRING=${all_file}
#       -D headers:STRING="${HEADERS}"
#       -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/GenerateAllHeaders.cmake
#  DEPENDS ${HEADERS}
#  COMMENT "Generating Mesquite_all_headers.hpp" )
#
SET_SOURCE_FILES_PROPERTIES(${all_file} PROPERTIES GENERATED true)

file( WRITE ${all_file} "#ifndef MESQUITE_ALL_HEADERS_HPP\n#define MESQUITE_ALL_HEADERS_HPP\n" )
foreach( header ${hdr_mesquite} )
  string( REGEX MATCH "[^/]+\\.h(pp)?" file ${header} )
  file( APPEND ${all_file} "#include \"${file}\"\n")
endforeach()
file( APPEND ${all_file} "#endif\n" )


#
# Set up package-specific options
#   - Sourced from <mesquite>/src/CMakeLists.txt
#

# Read the package version number specified as the second argument
# to the AC_INIT macro in a GNU Autoconf configure.in file.
#
# Define the following variables:
# VERSION_STRING:  The second argument to AC_INIT
# MAJOR_VERSION:   For a version string of the form m.n.p..., m
# MINOR_VERSION:   For a version string of the form m.n.p..., n
# PATCH_VERSION:   For a version string of the form m.n.p..., p

FILE( READ "configure.ac" configure_IN )
STRING( REGEX REPLACE "^.*AC_INIT *\\([^,]+, *([^, )]+).*$" "\\1" VERSION_STRING "${configure_IN}" )
STRING( REGEX REPLACE "^([0-9]+)(\\..*)?$" "\\1" MAJOR_VERSION "${VERSION_STRING}" )
STRING( REGEX REPLACE "^[0-9]+\\.([0-9]+)(\\..*)?$" "\\1" MINOR_VERSION "${VERSION_STRING}" )
IF ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )
  STRING( REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" PATCH_VERSION "${VERSION_STRING}" )
ELSE()
  SET(PATCH_VERSION "0")
ENDIF ( VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+" )

SET(MSQ_VERSION ${VERSION_STRING})
SET(MSQ_VERSION_MAJOR ${MAJOR_VERSION})
SET(MSQ_VERSION_MINOR ${MINOR_VERSION})
SET(MSQ_VERSION_PATCH ${PATCH_VERSION})
SET(MSQ_VERSION_STRING "\"Mesquite ${MSQ_VERSION}\"")
SET(Mesquite_NAMESPACE "Mesquite${MAJOR_VERSION}" CACHE STRING "Mesquite namespace" )
SET(MESQUITE_NS ${Mesquite_NAMESPACE})
IF (NOT MESQUITE_NS STREQUAL "Mesquite" )
  SET(MESQUITE_NS_ALIAS 1)
ENDIF()
MESSAGE(STATUS "Mesquite version: ${MSQ_VERSION}")

#check for specific functions
CHECK_FUNCTION_EXISTS(times     HAVE_TIMES)
CHECK_FUNCTION_EXISTS(clock     HAVE_CLOCK)
CHECK_FUNCTION_EXISTS(cbrt      HAVE_CBRT)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_FUNCTION_EXISTS(vsprintf  HAVE_VSPRINTF)
CHECK_FUNCTION_EXISTS(finite    HAVE_FINITE)
CHECK_FUNCTION_EXISTS(_Isfinite HAVE__ISFINITE)

#check for header files
CHECK_INCLUDE_FILE(ieeefp.h    MSQ_HAVE_IEEEFP_H)
CHECK_INCLUDE_FILE(ieeefp.h    MSQ_HAVE_STDINT_H)
CHECK_INCLUDE_FILE(inttypes.h  MSQ_HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(termios.h   HAVE_TERMIOS_H)
CHECK_INCLUDE_FILE(sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILE(dlfcn.h     HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(fcntl.h     HAVE_FCNTL_H)

IF(NOT HAVE_TIMES AND NOT HAVE_CLOCK)
  MESSAGE( FATAL_ERROR
    "\n***"
    "\n*** ERROR: Neither the function times nor the function clock could be found. At least one of these needs to exist or Mesquite cannot build properly."
    "\n***\n"
    )
ENDIF()

CONFIGURE_FILE(
  src/include/mesquite_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h
)
CONFIGURE_FILE(
  src/include/mesquite_version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_version.h
)


#
# Output config data
#   - Sourced from <mesquite>/src/CMakeLists.txt
#

MACRO(MESQUITE_APPEND_DEFINE DEFINE_LIST DEFINE_NAME)
  IF(${DEFINE_NAME})
    SET( ${DEFINE_LIST} "${${DEFINE_LIST}} -D${DEFINE_NAME}" )
  ENDIF()
ENDMACRO()

SET( MESQUITE_DEFINES "" )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES MSQ_TRAP_FPE )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES MSQ_ENABLE_DEBUG )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES NDEBUG )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES MSQ_USE_FUNCTION_TIMERS )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_TIMES )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_CLOCK )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_CBRT )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_VSNPRINTF )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_VSPRINTF )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_FINITE )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE__ISFINITE )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_IEEEFP_H )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_INTTYPES_H )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_TERMIOS_H )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_SYS_IOCTL_H )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_DLFCN_H )
MESQUITE_APPEND_DEFINE( MESQUITE_DEFINES HAVE_FCNTL_H )

SET_SOURCE_FILES_PROPERTIES( ${src_mesquite} COMPILE_FLAGS ${MESQUITE_DEFINES} )

#############################
# Alternative implementation?
#############################

# MACRO(MESQUITE_ADD_COMPILE_DEFINITION DEFINE_NAME)
#   IF(${DEFINE_NAME})
#     ADD_COMPILE_DEFINITIONS( "${DEFINE_NAME}" )
#   ENDIF()
# ENDMACRO()

# MESQUITE_ADD_COMPILE_DEFINITION( MSQ_TRAP_FPE )
# MESQUITE_ADD_COMPILE_DEFINITION( MSQ_ENABLE_DEBUG )
# MESQUITE_ADD_COMPILE_DEFINITION( NDEBUG )
# MESQUITE_ADD_COMPILE_DEFINITION( MSQ_USE_FUNCTION_TIMERS )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_TIMES )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_CLOCK )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_CBRT )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_VSNPRINTF )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_VSPRINTF )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_FINITE )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE__ISFINITE )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_IEEEFP_H )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_INTTYPES_H )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_TERMIOS_H )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_SYS_IOCTL_H )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_DLFCN_H )
# MESQUITE_ADD_COMPILE_DEFINITION( HAVE_FCNTL_H )

#
# Configure the library
#

DEAL_II_ADD_LIBRARY(obj_mesquite OBJECT ${src_mesquite} ${hdr_mesquite})