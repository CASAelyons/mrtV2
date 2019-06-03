SOURCES.c= run_mrtV2_alt.c mrtV2_linestring.c threshold_functions.c
$(CC)=gcc
CFLAGS=-Wall -I. -I/usr/local/include -I/usr/include/libxml2 -I/usr/local/libconfig/include -I/usr/local/jansson/include
LDLIBS= -L/usr/local/libconfig/lib -L/usr/local/jansson/lib -lm  -lnetcdf -ljansson -lconfig -lxml2

PROGRAM=mrtV2_linestring
OBJECTS= $(SOURCES.c:.c=.o)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:
	rm -f $(PROGRAM) $(OBJECTS)

install: all
	@echo "nothing to install"
