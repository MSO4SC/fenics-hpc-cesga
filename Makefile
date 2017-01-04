CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin` -g

DEST    = demo


OBJECTS = main.o MeshBC.o  LaplacianSmoother.o NodeNormal.o SpaceTimeFunction.o SlipBC.o  ufc2/NSEMomentum3D.o ufc2/NSEDensity3D.o ufc2/NSEContinuity3D.o ufc2/ProjectDensity3D.o ufc2/NSEResidualSC3D.o ufc2/NSEDualMomentum3D.o ufc2/NSEDualContinuity3D.o ufc2/NSEErrRepMomentum3D.o ufc2/NSEErrRepContinuity3D.o ufc2/Drag3D.o  ufc2/Laplacian2D.o  ufc2/Laplacian3D.o  ufc2/Poisson.o ufc2/L2ProjPfromM.o    ufc2/L2ProjUfromM2D.o  ufc2/L2ProjUfromM.o   



all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@
