DEST    = cube

all: $(DEST)

# this is the path to unicorn source
# it will be the path to install directory when install will be available
export LIBDIR=/opt/fenics-hpc/source/unicorn

include $(LIBDIR)/Makefile.lib

MYOBJECTS = cube.o

OBJECTS = $(UFC2OBJECTS) $(LIBOBJECTS) $(MYOBJECTS)

clean: cleanLib
	-rm -f $(MYOBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)
