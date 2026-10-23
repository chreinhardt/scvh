# Object files from the SCVH EOS library
SCVHEOS_OBJ = scvheos.o

OBJ = $(SCVHEOS_OBJ)

# The same objects built with -DSCVH_NOEXIT: failed lookups return
# SCVHEOS_FAIL_VALUE instead of exiting. Only the test suite uses them.
NOEXIT_OBJ = $(OBJ:.o=.noexit.o)

# Executables
EXE = testscvheos

# Definitions
DEFS =

# GNU Science library (uncomment if not needed)
GSL_LIB = -lgsl -lgslcblas

CFLAGS ?= -O3 $(DEFS) -Wall -std=c99

LIBS ?= -lm $(GSL_LIB)

default:
	@echo "Please specify which tool you want to make."

all: default

%.noexit.o: %.c
	$(CC) $(CFLAGS) -DSCVH_NOEXIT -c $< -o $@

#
# Test the SCVHEOS library.
#
# testscvheos is the test suite (invariants on any table, see testscvheos.c).
# It needs the non-exiting library for the failure-contract test, so it links
# scvheos.noexit.o, never the production scvheos.o. Run: ./testscvheos [iMat]
# with the table in the working directory; the exit status is the number of
# failed checks.
#
testscvheos: testscvheos.noexit.o $(NOEXIT_OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

check: testscvheos
	./testscvheos

testscvheosreadtable: testscvheosreadtable.o $(OBJ)
	cc -o testscvheosreadtable testscvheosreadtable.o $(OBJ) $(LIBS)

testscvheosinterp: testscvheosinterp.o $(OBJ)
	cc -o testscvheosinterp testscvheosinterp.o $(OBJ) $(LIBS)

testscvheospressure: testscvheospressure.o $(OBJ)
	cc -o testscvheospressure testscvheospressure.o $(OBJ) $(LIBS)

testscvheosintenergy: testscvheosintenergy.o $(OBJ)
	cc -o testscvheosintenergy testscvheosintenergy.o $(OBJ) $(LIBS)

testscvheosentropy: testscvheosentropy.o $(OBJ)
	cc -o testscvheosentropy testscvheosentropy.o $(OBJ) $(LIBS)

testscvheoscv: testscvheoscv.o $(OBJ)
	cc -o testscvheoscv testscvheoscv.o $(OBJ) $(LIBS)

testscvheosderivs: testscvheosderivs.o $(OBJ)
	cc -o testscvheosderivs testscvheosderivs.o $(OBJ) $(LIBS)

testscvheosinv: testscvheosinv.o $(OBJ)
	cc -o testscvheosinv testscvheosinv.o $(OBJ) $(LIBS)

testscvheossoundspeed: testscvheossoundspeed.o $(OBJ)
	cc -o testscvheossoundspeed testscvheossoundspeed.o $(OBJ) $(LIBS)

testscvheosunits: testscvheosunits.o $(OBJ)
	cc -o testscvheosunits testscvheosunits.o $(OBJ) $(LIBS)

scvheoscalcentropyonreos3grid: scvheoscalcentropyonreos3grid.o $(OBJ)
	cc -o scvheoscalcentropyonreos3grid scvheoscalcentropyonreos3grid.o $(OBJ) $(LIBS)

scvheos_h_limitedtoreos3: scvheos_h_limitedtoreos3.o $(OBJ)
	cc -o scvheos_h_limitedtoreos3 scvheos_h_limitedtoreos3.o $(OBJ) $(LIBS)

scvheos_he_limitedtoreos3: scvheos_he_limitedtoreos3.o $(OBJ)
	cc -o scvheos_he_limitedtoreos3 scvheos_he_limitedtoreos3.o $(OBJ) $(LIBS)

scvheos_hhe_limitedtoreos3: scvheos_hhe_limitedtoreos3.o $(OBJ)
	cc -o scvheos_hhe_limitedtoreos3 scvheos_hhe_limitedtoreos3.o $(OBJ) $(LIBS)

scvheoscalcentropy_gsl: scvheoscalcentropy_gsl.o $(OBJ)
	cc -o scvheoscalcentropy_gsl scvheoscalcentropy_gsl.o $(OBJ) $(LIBS)

test_ravit_model: test_ravit_model.o $(OBJ)
	cc -o test_ravit_model test_ravit_model.o $(OBJ) $(LIBS)

test_ravit_model_prhou: test_ravit_model_prhou.o $(OBJ)
	cc -o test_ravit_model_prhou test_ravit_model_prhou.o $(OBJ) $(LIBS)

extrap_scvh_on_grid: extrap_scvh_on_grid.o $(OBJ)
	cc -o extrap_scvh_on_grid extrap_scvh_on_grid.o $(OBJ) $(LIBS)

#
# Tools.
#
scvheoscalcentropy: scvheoscalcentropy.o $(OBJ)
	cc -o scvheoscalcentropy scvheoscalcentropy.o $(OBJ) $(LIBS)

scvheoscalcintenergy: scvheoscalcintenergy.o $(OBJ)
	cc -o scvheoscalcintenergy scvheoscalcintenergy.o $(OBJ) $(LIBS)

scvheosprintmat: scvheosprintmat.o $(OBJ)
	cc -o scvheosprintmat scvheosprintmat.o $(OBJ) $(LIBS)

scvheoscalctemprhou: scvheoscalctemprhou.o $(OBJ)
	cc -o scvheoscalctemprhou scvheoscalctemprhou.o $(OBJ) $(LIBS)

scvheoscalcurhotemp: scvheoscalcurhotemp.o $(OBJ)
	cc -o scvheoscalcurhotemp scvheoscalcurhotemp.o $(OBJ) $(LIBS)

scvheos_calc_pofrhot: scvheos_calc_pofrhot.o $(OBJ)
	cc -o scvheos_calc_pofrhot scvheos_calc_pofrhot.o $(OBJ) $(LIBS)

scvheos_calc_rhoofpt: scvheos_calc_rhoofpt.o $(OBJ)
	cc -o scvheos_calc_rhoofpt scvheos_calc_rhoofpt.o $(OBJ) $(LIBS)

calc_model_entropy: calc_model_entropy.o $(OBJ)
	cc -o calc_model_entropy calc_model_entropy.o $(OBJ) $(LIBS)

calc_ravit_model_entropy: calc_ravit_model_entropy.o $(OBJ)
	cc -o calc_ravit_model_entropy calc_ravit_model_entropy.o $(OBJ) $(LIBS)

calc_ravit_model_codeunits: calc_ravit_model_codeunits.o $(OBJ)
	cc -o calc_ravit_model_codeunits calc_ravit_model_codeunits.o $(OBJ) $(LIBS)

scvheos_calc_isentrope: scvheos_calc_isentrope.o $(OBJ)
	cc -o scvheos_calc_isentrope scvheos_calc_isentrope.o $(OBJ) $(LIBS)

calc_isentrope_rhot: calc_isentrope_rhot.o $(OBJ)
	cc -o calc_isentrope_rhot calc_isentrope_rhot.o $(OBJ) $(LIBS)

fix_lowt_entropy: fix_lowt_entropy.o $(OBJ)
	cc -o fix_lowt_entropy fix_lowt_entropy.o $(OBJ) $(LIBS)

calc_rhot_isentropes: calc_rhot_isentropes.o $(OBJ)
	cc -o calc_rhot_isentropes calc_rhot_isentropes.o $(OBJ) $(LIBS)

.PHONY: default all check clean cleanall

clean:
	rm -f $(OBJ) $(NOEXIT_OBJ) testscvheos.noexit.o

cleanall:
	rm -f $(EXE) $(OBJ) $(NOEXIT_OBJ) *.o