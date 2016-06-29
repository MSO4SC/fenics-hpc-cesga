FFCFLAGS = -O -r quadrature -fsplit -l dolfin

all:
	ffc $(FFCFLAGS) NSEMomentum3D.ufl
	ffc $(FFCFLAGS) NSEContinuity3D.ufl
	ffc $(FFCFLAGS) NSEDualMomentum3D.ufl 
	ffc $(FFCFLAGS) NSEDualContinuity3D.ufl 
	ffc $(FFCFLAGS) NSEErrRepMomentum3D.ufl 
	ffc $(FFCFLAGS) NSEErrRepContinuity3D.ufl 
	ffc $(FFCFLAGS) Drag3D.ufl 
	ffc $(FFCFLAGS) NSEH1.ufl 
	ffc $(FFCFLAGS) NSEH12.ufl 
	ffc $(FFCFLAGS) NSEMomentumResidual3D.ufl
	ffc $(FFCFLAGS) NSEContinuityResidual3D.ufl
	ffc $(FFCFLAGS) NSEH1Momentum3D.ufl
	ffc $(FFCFLAGS) NSEH1Continuity3D.ufl
	ffc $(FFCFLAGS) NSEH1MomentumGlobal3D.ufl
	ffc $(FFCFLAGS) NSEH1ContinuityGlobal3D.ufl
	ffc $(FFCFLAGS) NSEMomentumResidualGlobal3D.ufl
	ffc $(FFCFLAGS) NSEContinuityResidualGlobal3D.ufl
	ffc $(FFCFLAGS) NSEErrEst.ufl
	ffc $(FFCFLAGS) NSEErrEstGlobal.ufl
