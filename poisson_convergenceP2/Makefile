CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = poisson
OBJECTS = main.o DirichletBC_CG2v_CG2.o

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
