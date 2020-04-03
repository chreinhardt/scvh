
# Object files from the SCVH EOS library
SCVHEOS_OBJ = scvheos.o

OBJ = $(SCVHEOS_OBJ)

# Executables
EXE =

# Definitions
DEFS =

# GNU Science library (uncomment if not needed)
#GSL_LIB = -lgsl -lgslcblas

CFLAGS ?= -O3 $(defs) -Wall -std=c99

LIBS ?= -lm $(GSL_LIB)

default:
	@echo "Please specify which tool you want to make."

all:
	default

#
# Test the SCVHEOS library.
#
testscvheos: testscvheos.o $(OBJ)
	cc -o testscvheos testscvheos.o $(OBJ) $(LIBS)

testscvheosreadtable: testscvheosreadtable.o $(OBJ)
	cc -o testscvheosreadtable testscvheosreadtable.o $(OBJ) $(LIBS)


clean:
	rm $(OBJ)

cleanall:
	rm $(EXE) $(OBJ) *.o
