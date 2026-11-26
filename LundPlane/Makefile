# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=`$(FASTJETCONFIG) --prefix`
CXX=g++
CXXFLAGS= -O3 -Wall -g
install_script = $(SHELL) ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=LundPlane
SRCS=LundGenerator.cc LundWithSecondary.cc SecondaryLund.cc RecursiveLundEEGenerator.cc
EXAMPLES=example example_secondary example_dpsi_collinear example_dpsi_slice
INSTALLED_HEADERS=LundGenerator.hh LundWithSecondary.hh SecondaryLund.hh RecursiveLundEEGenerator.hh LundJSON.hh LundEEHelpers.hh LundPlane.hh
# List here any dependencies on other contribs. Note that they 
# must be listed (potentially with versioning constraints) also
# in the FJCONTRIB.cfg file.
DEPENDS_ON =
#------------------------------------------------------------------------

CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
LDFLAGS += -lm $(shell $(FASTJETCONFIG) --libs --plugins)

INCLUDE += -Iinclude/
# loop over all entries in DEPENDS_ON and for each entry ENTRY
# add -I../ENTRY to INCLUDE, -L../ENTRY to LDFLAGS and -lENTRY to LDFLAGS
INCLUDE += $(foreach ENTRY,$(DEPENDS_ON),-I../$(ENTRY)/include)
LDFLAGS += $(foreach ENTRY,$(DEPENDS_ON),-L../$(ENTRY) -l$(ENTRY))
CXXFLAGS+= $(INCLUDE)

OBJS  = $(SRCS:.cc=.o)
EXAMPLES_SRCS  = $(EXAMPLES:=.cc)

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d
install_DATA    = $(install_script) -c -m 644
install_PROGRAM = $(install_script) -c -s
install_SCRIPT  = $(install_script) -c

.PHONY: clean distclean examples check install

# compilation of the code (default target)
all: lib$(NAME).a

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

# building the examples
examples: $(EXAMPLES)

# the following construct makes it possible to automatically build
# each of the examples listed in $EXAMPLES
$(EXAMPLES): % : %.o all
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -l$(NAME) $(LDFLAGS)

# check that everything went fine
check: examples
	@for prog in $(EXAMPLES); do\
	  if [ "$${prog}" = "example_dpsi_collinear" ] || [ "$${prog}" = "example_dpsi_slice" ]; then \
	    $(check_script) $${prog} ../data/single-ee-event.dat || exit 1; \
	  else \
	    $(check_script) $${prog} ../data/single-event.dat || exit 1; \
	  fi; \
	done
	@echo "All tests successful"

# cleaning the directory
clean:
	rm -f *~ *.o

distclean: clean
	rm -f lib$(NAME).a $(EXAMPLES)

# install things in PREFIX/...
install: all
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	for header in $(INSTALLED_HEADERS); do\
	  $(install_HEADER) include/fastjet/contrib/$$header $(PREFIX)/include/fastjet/contrib/;\
	done
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).a $(PREFIX)/lib

depend:
	makedepend -Y -- -Iinclude -- $(SRCS) $(EXAMPLES_SRCS)
# DO NOT DELETE

LundGenerator.o: include/fastjet/contrib/LundGenerator.hh
LundWithSecondary.o: include/fastjet/contrib/LundWithSecondary.hh
LundWithSecondary.o: include/fastjet/contrib/LundGenerator.hh
LundWithSecondary.o: include/fastjet/contrib/SecondaryLund.hh
SecondaryLund.o: include/fastjet/contrib/SecondaryLund.hh
SecondaryLund.o: include/fastjet/contrib/LundGenerator.hh
RecursiveLundEEGenerator.o: include/fastjet/contrib/RecursiveLundEEGenerator.hh
RecursiveLundEEGenerator.o: include/fastjet/contrib/LundEEHelpers.hh
example.o: include/fastjet/contrib/LundGenerator.hh
example.o: include/fastjet/contrib/LundJSON.hh
example.o: include/fastjet/contrib/LundGenerator.hh
example_secondary.o: include/fastjet/contrib/LundWithSecondary.hh
example_secondary.o: include/fastjet/contrib/LundGenerator.hh
example_secondary.o: include/fastjet/contrib/SecondaryLund.hh
example_dpsi_collinear.o: include/fastjet/contrib/RecursiveLundEEGenerator.hh
example_dpsi_collinear.o: include/fastjet/contrib/LundEEHelpers.hh
example_dpsi_slice.o: include/fastjet/contrib/RecursiveLundEEGenerator.hh
example_dpsi_slice.o: include/fastjet/contrib/LundEEHelpers.hh
