CFLAGS  = `pkg-config --cflags dolfin` -I./ -I./ufc2
LIBS    = `pkg-config --libs dolfin`  -g
CXX     = `pkg-config --variable=compiler dolfin`

TARGET  = cube
FORM_SRC	 = compressible2d_p.cpp L2proj2D_momen.cpp L2proj2D_energ.cpp L2proj2D_densi.cpp compressible2d_mach.cpp compressible2d_sound.cpp compressible2d_velo.cpp compressible2d.cpp DirichletBC_CG1v_CG1.cpp L2proj2D.cpp

FORM_OBJ = ${FORM_SRC:.cpp=.o}

all: $(TARGET)

install:

clean:
	-rm -f *.o core *.core main.o $(FORM_OBJ)  $(TARGET)

forms:	$(FORM_OBJ)

solver: $(SOLVER_SRC)
	$(MAKE) -C ../

$(TARGET): forms
	$(CXX) $(CFLAGS) main.cpp $(FORM_OBJ)  $(LIBS) -o $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

