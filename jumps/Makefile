CFLAGS  = `pkg-config --cflags dolfin` -I./ -I../
LIBS    = `pkg-config --libs dolfin`  -g
CXX     = `pkg-config --variable=compiler dolfin`

TARGET  = cube
FORM_SRC	 = basic2D.cpp

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

