FFCFLAGS = -O -r quadrature -fsplit -l dolfin

all:
	ffc $(FFCFLAGS) Drag3D.ufl
	ffc $(FFCFLAGS) Laplacian2D.ufl
	ffc $(FFCFLAGS) Laplacian3D.ufl
	ffc $(FFCFLAGS) NSEbase.ufl
	ffc $(FFCFLAGS) NSEContinuity2D.ufl
	ffc $(FFCFLAGS) NSEContinuity3D.ufl
	ffc $(FFCFLAGS) NSEDensity3D.ufl
	ffc $(FFCFLAGS) NSEDualContinuity3D.ufl
	ffc $(FFCFLAGS) NSEDualMomentum3D.ufl
	ffc $(FFCFLAGS) NSEErrRepContinuity3D.ufl
	ffc $(FFCFLAGS) NSEErrRepMomentum3D.ufl
	ffc $(FFCFLAGS) NSEMomentum2D.ufl
	ffc $(FFCFLAGS) NSEMomentum3D.ufl
	ffc $(FFCFLAGS) NSEResidualSC3D.ufl
	ffc $(FFCFLAGS) Poisson.ufl
	ffc $(FFCFLAGS) ProjectDensity3D.ufl
	ffc $(FFCFLAGS) L2ProjUfromM2D.ufl  
	ffc $(FFCFLAGS) L2ProjPfromM.ufl    
	ffc $(FFCFLAGS) L2ProjUfromM.ufl   