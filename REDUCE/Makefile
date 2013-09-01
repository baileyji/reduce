#
#  Platform independant makefile for the SME Call_External's.
#
SHELL=/bin/sh
CFLAGS=
FFLAGS=
F_LD_POST=

# The following is the default entry point. This section will determine 
# what system we are on, set the correct flags and call this same makefile
# again with the correct flags.

all : 
	@echo "OS type detected: "`uname`
	@case `uname` in \
           "SunOS") if [  -d /proc ]; then \
                        make libs \
                        "SHL_EXT=so.sunos.sparc.32" \
                        "CC=CC" \
                        "C_LD=CC" \
                        "CFLAGS=-Kpic -DSPARC -G -O3" \
                        "FFLAGS=-pic -G -O3 -ext_names=plain" \
                        "C_LD_FLAGS=-G -DSOLARIS" \
                        "C_LD_POST= -lsunmath -lF77 -lM77 -lm -lc"; \
                    fi  \
                ;; \
           "SunOS64") if [  -d /proc ]; then \
                        make libs \
                        "SHL_EXT=so.sunos.sparc.64" \
                        "CC=CC" \
                        "C_LD=CC" \
                        "CFLAGS=-xtarget=ultra -xarch=v9 -Kpic -G -O3" \
                        "C_LD_FLAGS=-xtarget=ultra -xarch=v9 -G" \
                        "C_LD_POST= -lsunmath -lF77 -lM77 -lm -lc"; \
                    fi  \
                ;; \
	   "HP-UX") make libs \
			"SHL_EXT=so.hp-ux.risc.32" \
			"CC=aCC" \
			"C_LD=ld" \
			"CFLAGS=+z +O2 +Onolimit +Oloop_unroll" \
			"C_LD_FLAGS=-b" \
			"C_LD_POST= -lm -lc";; \
	   "Linux32" ) make libs \
			"SHL_EXT=so.linux.x86_32.32" \
			"CC=cc" \
			"F77=f77 -O" \
			"C_LD=cc" \
			"F_LD=f77" \
			"CFLAGS=-m32 -shared -O2 -fPIC" \
			"C_LD_FLAGS=-m32 -O2 -shared -fPIC" \
			"C_LD_POST= ";; \
	   "Linux" ) make libs \
			"SHL_EXT=so.linux.x86_64.64" \
			"CC=cc" \
			"F77=f77 -O" \
			"C_LD=cc" \
			"F_LD=f77" \
			"CFLAGS=-m64 -shared -O2 -fPIC" \
			"C_LD_FLAGS=-m64 -O2 -shared -fPIC" \
			"C_LD_POST=";; \
	   "Darwin 32" ) make libs \
			"SHL_EXT=so."$(OSTYPE)"."$(MACHTYPE)".32" \
			"CC=cc" \
			"C_LD=cc" \
			"CFLAGS=-x c++ -fast -fPIC -no-cpp-precomp -dynamic -fPIC -fno-common -D_REENTRANT" \
			"C_LD_FLAGS=-fast -bundle -flat_namespace -undefined suppress" \
			"C_LD_POST= ";; \
           "Darwin" ) make libs \
                        "SHL_EXT=so."$(OSTYPE)"."$(MACHTYPE)".64" \
                        "CC=cc" \
                        "C_LD=cc" \
                        "CFLAGS=-m64 -fast -dynamic -fPIC" \
                        "C_LD_FLAGS= -m64 -fast -dynamiclib -i-static " \
                        "C_LD_POST= ";; \
 	   *) echo "This system is not supported" ;; \
       esac

libs :   bandsol.$(SHL_EXT)

# C Only libs

bandsol.$(SHL_EXT) : bandsol.o Makefile
	$(C_LD) $(C_LD_FLAGS) -o bandsol.$(SHL_EXT) bandsol.o $(C_LD_POST)

bandsol.o: bandsol.c Makefile
	$(CC) $(CFLAGS) -c bandsol.c
