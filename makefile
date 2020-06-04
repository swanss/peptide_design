MSTDIR = ../MST
MSTINCL = $(MSTDIR)/include
MSTSRC = $(MSTDIR)/src
MSTOBJS = $(MSTDIR)/objs
MSTLIB = $(MSTDIR)/lib

STRUCTGENDIR = ../structgen
STRUCTGENINCL = $(STRUCTGENDIR)/include
STRUCTGENSRC = $(STRUCTGENDIR)/src
STRUCTGENOBJS = $(STRUCTGENDIR)/objs
STRUCTGENLIB = $(STRUCTGENDIR)/lib

CC = g++
CFLAGS = -std=c++11 -g -gdwarf-3 -O3 -fPIC -I$(MSTINCL) -I$(STRUCTGENINCL) -I$(INCL) # -g -gdwarf-3
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

LIB         := -L$(MSTLIB) -L$(STRUCTGENLIB) -lmst -lstructgen -lmstcondeg -lmstfuser -lmstoptim -lmstfasst -lmstmagic -ldtermen
INC         :=
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
PYLIB_PATH = $(shell python3-config --exec-prefix)/lib
PYLIB = -L$(PYLIB_PATH) $(LIB) -L$(OUT)/lib -lsandbox -L/home/ifs-users/venkats/miniconda3/include/python3.7m -L/usr/local/lib -L/home/ifs-users/venkats/lib -lboost_python37 $(shell python3-config --libs) -Wl,--no-undefined
PYFLAGS = $(shell python3-config --includes) -I/home/ifs-users/venkats/include -O2 -fPIC -std=c++11 -I$(INCL) -I$(MSTINCL) -I$(STRUCTGENINCL) -L$(OUT)/lib

# make the boost.python shared object
python: $(OUT)/lib/sandbox.so

$(BUILDDIR)/python.o: $(SRCDIR)/python.cpp # $(OBJECTS)
	$(CC) $(PYFLAGS) -c -o $@ $<

$(OUT)/lib/sandbox.so: $(BUILDDIR)/python.o $(OBJECTS) lib
	$(CC) -shared -o $@ $< $(PYLIB) -Wl,-rpath,$(PYLIB_PATH)

.PHONY: lib
lib: $(OUT)/lib/libsandbox.a

$(OUT)/lib/libsandbox.a: directories $(OBJECTS)
	# Remove the Python.o file (it will need to be rebuilt anyway)
	# rm -f $(BUILDDIR)/python.o

	ar rs $(OUT)/lib/libsandbox.a $(OBJECTS) $(MSTOBJS)/*.o $(STRUCTGENOBJS)/*.o

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
	$(CC) $(CFLAGS) $^ $(LIB) -o $@

$(BIN)/%: $(TEST)/%.$(SRCEXT) $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LIB) -o $@

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT) directories
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LIB) -c -o $@ $<
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
