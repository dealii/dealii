SRC_MISC_DIR = Misc
MSQ_INCLUDES += -I$(top_srcdir)/src/$(SRC_MISC_DIR)
INCLUDED_MAKEFILES += $(SRC_MISC_DIR)/Makefile.am

if ENABLE_ITAPS
  MSQ_ITAPS_COMMON_SRC = $(SRC_MISC_DIR)/MsqIBase.cpp \
                         $(SRC_MISC_DIR)/MsqIBase.hpp
else
  MSQ_ITAPS_COMMON_SRC =
endif

MSQ_SRCS += \
	$(SRC_MISC_DIR)/BoundedCylinderDomain.cpp \
	$(SRC_MISC_DIR)/ConicDomain.cpp \
	$(SRC_MISC_DIR)/CylinderDomain.cpp \
	$(SRC_MISC_DIR)/DomainClassifier.cpp \
	$(SRC_MISC_DIR)/DomainUtil.cpp \
	$(SRC_MISC_DIR)/DomainUtil.hpp \
	$(SRC_MISC_DIR)/EdgeIterator.cpp \
	$(SRC_MISC_DIR)/Exponent.cpp \
	$(SRC_MISC_DIR)/FileTokenizer.cpp \
	$(SRC_MISC_DIR)/MeshDomain1D.cpp \
	$(SRC_MISC_DIR)/MeshTransform.cpp \
	$(SRC_MISC_DIR)/MesquiteVersion.cpp \
	$(SRC_MISC_DIR)/MsqDebug.cpp \
	$(SRC_MISC_DIR)/MsqError.cpp \
	$(SRC_MISC_DIR)/MsqFPE.cpp \
	$(SRC_MISC_DIR)/MsqGeomPrim.cpp \
	$(SRC_MISC_DIR)/MsqHessian.cpp \
	$(SRC_MISC_DIR)/MsqInterrupt.cpp \
	$(SRC_MISC_DIR)/MsqTimer.cpp \
	$(SRC_MISC_DIR)/PlanarDomain.cpp \
	$(SRC_MISC_DIR)/SimpleStats.cpp \
	$(SRC_MISC_DIR)/SphericalDomain.cpp \
	$(SRC_MISC_DIR)/Vector3D.cpp \
	$(SRC_MISC_DIR)/VtkTypeInfo.cpp \
	$(SRC_MISC_DIR)/XYPlanarDomain.cpp \
	$(SRC_MISC_DIR)/XYRectangle.cpp \
        $(MSQ_ITAPS_COMMON_SRC)


MSQ_HDRS += \
	$(SRC_MISC_DIR)/BoundedCylinderDomain.hpp \
	$(SRC_MISC_DIR)/Bits.hpp \
	$(SRC_MISC_DIR)/ConicDomain.hpp \
	$(SRC_MISC_DIR)/CylinderDomain.hpp \
	$(SRC_MISC_DIR)/DomainClassifier.hpp \
	$(SRC_MISC_DIR)/EdgeIterator.hpp \
	$(SRC_MISC_DIR)/Exponent.hpp \
	$(SRC_MISC_DIR)/FileTokenizer.hpp \
	$(SRC_MISC_DIR)/Matrix3D.hpp \
	$(SRC_MISC_DIR)/MeshDomain1D.hpp \
	$(SRC_MISC_DIR)/MeshTransform.hpp \
	$(SRC_MISC_DIR)/MsqDebug.hpp \
	$(SRC_MISC_DIR)/MsqError.hpp \
	$(SRC_MISC_DIR)/MsqFPE.hpp \
	$(SRC_MISC_DIR)/MsqGeomPrim.hpp \
	$(SRC_MISC_DIR)/MsqHessian.hpp \
	$(SRC_MISC_DIR)/MsqInterrupt.hpp \
	$(SRC_MISC_DIR)/MsqMatrix.hpp \
	$(SRC_MISC_DIR)/MsqTimer.hpp \
	$(SRC_MISC_DIR)/PlanarDomain.hpp \
	$(SRC_MISC_DIR)/SimpleStats.hpp \
	$(SRC_MISC_DIR)/SphericalDomain.hpp \
	$(SRC_MISC_DIR)/SymMatrix3D.hpp \
	$(SRC_MISC_DIR)/Vector3D.hpp \
	$(SRC_MISC_DIR)/VtkTypeInfo.hpp \
	$(SRC_MISC_DIR)/XYPlanarDomain.hpp \
	$(SRC_MISC_DIR)/XYRectangle.hpp

