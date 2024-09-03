#  makefile

# -- directories for the source code, object files, executable files, module files, and library files
SRC  = ./src/
OBJ  = ./obj/
BIN  = ./bin/
MOD  = ./mod/
LIB  = ./lib/
DIRS =	 \
	$(SRC) \
	$(OBJ) \
	$(BIN) \
	$(MOD) \
	$(LIB) \

# -- source code suffix, assuming all the same
SUF = f90
#SUF = f08

# -- the program name and executable to generate
PROGRAM = EI_Exc_Calc
PROGBIN = $(BIN)$(PROGRAM).e

# -- lists of source files, object files, and mod files.
SRCS = $(wildcard $(SRC)*.$(SUF))
OBJS = $(patsubst $(SRC)%,$(OBJ)%,$(patsubst %.$(SUF),%.o,$(SRCS)))
MODS = $(patsubst $(MOD)$(PROGRAM).mod,,$(patsubst $(SRC)%,$(MOD)%,$(patsubst %.$(SUF),%.mod,$(SRCS))))

# -- compiler
FC      = gfortran

# -- warning flags
FWFLAGS = \
	-Werror \
	-Wno-tabs \
	-Wno-unused-function \
	-Wno-maybe-uninitialized \
	-Wno-conversion
	#-Wall 

# -- compiler flags
FCFLAGS = \
	-J$(MOD) \
	-g \
	-ffpe-trap=zero,overflow,invalid \
	-fcheck=all,bounds \
	-fimplicit-none \
	-ffree-line-length-0 \
	-floop-parallelize-all \
	-mcmodel=medium \
	-pg 

# -- linker flags
FLFLAGS = \
	-I$(MOD) \
	-g \
	-fbacktrace \
	-L$(LIB) \
  -lstdlib \
  -llapack \
	-lopenblas 
	#-pg

# -- optimization flags
FOFLAGS = -Og

.PHONY: default debug clean superclean $(DIRS)

# -- parameterless make defaults to first option
default: $(DIRS) $(PROGBIN)

# -- make directories
$(DIRS) :
	@mkdir -p $@

# -- compile code, generate object and mod files.
#    If a source file changes, recompile only its object file. If the Makefile change, recompile everything
$(OBJS) : $(OBJ)%.o : $(SRC)%.$(SUF) Makefile
	$(FC) $(FCFLAGS) $(FWFLAGS) $(FOFLAGS) -c $< -o $@

# -- linker
$(PROGBIN) : $(OBJS)
	$(FC) -o $(PROGBIN) $^ $(FLFLAGS)

# -- for when stuff doesn't work
debug:
	@echo ""
	@echo "---------- debug ----------"
	@echo ""
	@echo "DIRS = $(DIRS)"
	@echo "DIRS = $(DIRS)"
	@echo ""
	@echo "SRCS = $(SRCS)"
	@echo ""
	@echo "OBJS = $(OBJS)"
	@echo ""
	@echo "MODS = $(MODS)"
	@echo ""
	@echo "PROGRAM = $(PROGRAM)"
	@echo ""
	@echo "PROGBIN = $(PROGBIN)"
	@echo ""
	@echo "ldd $(PROGBIN):"
	@ldd  "$(PROGBIN)"
	@echo ""
	@echo "FCFLAGS = $(FCFLAGS)"
	@echo ""
	@echo "FLFLAGS = $(FLFLAGS)"
	@echo ""
	@echo "FOFLAGS = $(FOFLAGS)"
	@echo ""
	@echo "---------------------------"
	@echo ""

# -- clean up main directory and object directory
clean:
	rm -rf $(OBJ)*.o $(PROGBIN) *.mod *.o *.out fort*

# -- clean but also remove mod files from module directory
superclean: clean
	rm -rf $(MOD)*.mod

# -- end of makefile
