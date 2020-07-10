MSTDIR = ../MST
MSTINCL = $(MSTDIR)/include
MSTSRC = $(MSTDIR)/src
MSTOBJS = $(MSTDIR)/objs
MSTLIB = $(MSTDIR)/lib

CC = g++
CFLAGS = -std=c++11 -g -gdwarf-3 -O3 -fPIC -I$(MSTINCL) -I$(INCL) # -g -gdwarf-3
MPICC = mpic++
MPIFLAGS = -std=c++0x -O3

OUT = .
INCL = $(OUT)/inc
SRC = $(OUT)/src
LIB = $(OUT)/lib
OBJS = $(OUT)/objs
TEST = $(OUT)/tests
TESTFILES = $(OUT)/testfiles
SCRIPTS = $(OUT)/scripts
PROGS = $(OUT)/programs
BIN = $(OUT)/bin
SENTINEL = .dir_sentinel

SRCDIR      := $(OUT)/src
INCDIR      := $(OUT)/include
BUILDDIR    := $(OUT)/objs
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

LIBFLAGS    := -L$(MSTLIB) -lmst -lmstcondeg -lmstfuser -lmstoptim -lmstfasst -lmstmagic -ldtermen
INCDEP      :=

SOURCES     := $(shell find $(SRCDIR) -not -path '*/\.*' -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
PROGRAMS 	:= $(patsubst $(PROGS)/%.$(SRCEXT),$(BIN)/%,$(shell find $(PROGS) -not -path '*/\.*' -type f -name *.$(SRCEXT)))
TESTPROGS 	:= $(patsubst $(TEST)/%.$(SRCEXT),$(BIN)/%,$(shell find $(TEST) -not -path '*/\.*' -type f -name *.$(SRCEXT)))

# Filter out python
SOURCES     := $(filter-out $(SRCDIR)/python.cpp, $(SOURCES))
OBJECTS     := $(filter-out $(BUILDDIR)/python.o, $(OBJECTS))

##############
# Main Targets
##############

.PHONY: all
all: $(PROGRAMS) $(TESTPROGS)

.PHONY: test
test: $(TESTPROGS)

.PHONY: remake
remake: cleaner all

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
	rm -rf $(BIN)
	rm -rf $(OUT)/lib


######################
# Boost Python
######################

uname := $(shell uname -s)
PYTHON_SUFFIX = $(shell $(pythonExec) -c "import sys; print(''.join(map(str,sys.version_info[0:2])));")
ifeq ($(uname),Linux)
	# TODO verify that this works on Linux
	pythonExec := python3
	PYLIB_PATH = $(shell $(pythonExec)-config --exec-prefix)/lib64
	PYLIB = -L$(PYLIB_PATH) $(LIBFLAGS) -L$(LIB) -ldl -lboost_python$(PYTHON_SUFFIX) -lpeptide_design $(LDLIBS)
else # MacOS
	# Requires python 3.8
	pythonExec := python3.8
	PYLIB_PATH = $(shell $(pythonExec) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'));")
	PYLIB = -L$(PYLIB_PATH) $(LIBFLAGS) -L$(LIB) -ldl -framework CoreFoundation -undefined dynamic_lookup -lboost_python$(PYTHON_SUFFIX) -lpeptide_design $(LDLIBS)
endif
PY_INCLUDES = $(shell $(pythonExec)-config --includes)
PY_SITE_INCLUDE_PARENT = $(shell $(pythonExec)-config --exec-prefix)
PYFLAGS = $(PY_INCLUDES) -I$(PY_SITE_INCLUDE_PARENT)/include -O3 -fPIC -std=c++11 -I$(INCL) -I$(MSTINCL)


# make the boost.python shared object
python: $(OUT)/lib/peptide_design.so

$(BUILDDIR)/python.o: $(SRCDIR)/python.cpp # $(OBJECTS)
	$(CC) $(PYFLAGS) -c -o $@ $<

$(OUT)/lib/peptide_design.so: $(BUILDDIR)/python.o $(OBJECTS) lib
	$(CC) -shared -o $@ $< $(PYLIB) -Wl,-rpath,$(PYLIB_PATH)

.PHONY: lib
lib: $(OUT)/lib/libpeptide_design.a

$(OUT)/lib/libpeptide_design.a: directories $(OBJECTS)
	# Remove the Python.o file (it will need to be rebuilt anyway)
	# rm -f $(BUILDDIR)/python.o

	ar rs $(OUT)/lib/libpeptide_design.a $(OBJECTS) $(MSTOBJS)/*.o

######################
# Make the Directories
######################

directories: $(BIN)/$(SENTINEL) $(BUILDDIR)/$(SENTINEL) $(OUT)/lib/$(SENTINEL)
	touch $@

$(BIN)/$(SENTINEL):
	mkdir -p $(BIN)
	touch $(BIN)/$(SENTINEL)

$(BUILDDIR)/$(SENTINEL):
	mkdir -p $(BUILDDIR)
	touch $(BUILDDIR)/$(SENTINEL)

$(OUT)/lib/$(SENTINEL):
	mkdir -p $(OUT)/lib
	touch $(OUT)/lib/$(SENTINEL)

#######################
# Building Source Files
#######################


-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

$(BIN)/%: $(PROGS)/%.$(SRCEXT) $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LIBFLAGS) -o $@

$(BIN)/%: $(TEST)/%.$(SRCEXT) $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LIBFLAGS) -o $@

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT) directories
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp


#######
# Debug
#######

.PHONY: print-%
print-%  : ; @echo $* = $($*)
