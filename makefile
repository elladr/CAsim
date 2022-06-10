.SUFFIXES: .f90

# The files including modules should be listed first

SOURCES= m_dengrowth.f90     simParam.f90  \
         dengrowth.f90         main.f90   \
         output.f90         
            

EXEC= a.out

# compiler for Fortran
FC = gfortran

FFLAGS = -O3 -shared  
LDFLAGS =
LIBS = 


OBJECTS = $(SOURCES:.f90=.o) $(LIBS)


$(EXEC): $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(LDFLAGS)


#This file is the first file in the SOURCES files 
m_dengrowth.o: m_dengrowth.f90
	gfortran -c m_dengrowth.f90

.f90.o: 
	$(FC) -c $<
