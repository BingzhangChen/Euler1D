.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c $(FFLAGS) $*.f90 -o $*.o
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

MAIN = Euler1D
#Source file
SRCS = Readcsv.f90\
    tridiagonal.f90\
	params.f90\
	Trait_functions.f90\
        Calc_PAR.f90\
	grid.f90\
	gridinterp.f90\
	time_interp.f90\
	variables.f90\
	Time_settings.f90\
	forcing.f90\
	netcdf_IO.f90\
 	Advection_center.f90\
	Diff_center.f90\
	GMK98_Size.f90\
	GMK98_simple.f90\
	Geider_Euler.f90\
	initialization.f90\
	timestep.f90\
	Main.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))

#Executable file
$(MAIN): $(OBJS)
	$(CFT) $(FFLAGS) $(LFLAGS) -o $(MAIN) $(OBJS)
