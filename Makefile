DEST    = demo

all: $(DEST)

# this is the path to unicorn source
# it will be the path to install directory when install will be available
export LIBDIR=/opt/fenics-hpc/source/unicorn

include $(LIBDIR)/Makefile.lib

MYOBJECTS = main.o

OBJECTS = $(UFC2OBJECTS) $(LIBOBJECTS)
$(MYOBJECTS) : $(OBJECTS)

clean: cleanLib
	-rm -f $(MYOBJECTS) $(DEST)

$(DEST): $(MYOBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(MYOBJECTS) $(CFLAGS) $(LIBS)
