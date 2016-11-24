CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin` -g

DEST    = cube
OBJECTS = cube.o ../NodeNormal.o ../SpaceTimeFunction.o ../SlipBC.o ../ufc2/NSEMomentum3D.o ../ufc2/NSEContinuity3D.o ../ufc2/NSEDualMomentum3D.o ../ufc2/NSEDualContinuity3D.o ../ufc2/NSEErrRepMomentum3D.o ../ufc2/NSEErrRepContinuity3D.o ../ufc2/Drag3D.o ../ufc2/NSEH1.o ../ufc2/NSEH12.o ../ufc2/NSEErrEst.o ../ufc2/NSEErrEstGlobal.o ../ufc2/NSEMomentumResidual3D.o ../ufc2/NSEContinuityResidual3D.o ../ufc2/NSEMomentumResidualGlobal3D.o ../ufc2/NSEContinuityResidualGlobal3D.o ../ufc2/NSEH1Momentum3D.o ../ufc2/NSEH1Continuity3D.o ../ufc2/NSEH1MomentumGlobal3D.o ../ufc2/NSEH1ContinuityGlobal3D.o

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@
