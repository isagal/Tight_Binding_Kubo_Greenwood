CXX = c++

# C++ flags
CFLAGS = -O3 -std=c++17
MKLROOT  = /opt/intel/oneapi/mkl
INCLUDES = -I $(MKLROOT)/2023.0.0/include/ 
LIBS     =  -L$(MKLROOT)/2023.0.0/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# executable file
EXECUTABLE = run

all: run
		
main.o: main.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c main.cpp
IO_data.o: IO_data.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c IO_data.cpp		
Watch.o: Watch.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Watch.cpp		
Hamiltonian.o: Hamiltonian.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Hamiltonian.cpp	
Psy.o: Psy.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Psy.cpp	
KuboGreenwood.o: KuboGreenwood.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c KuboGreenwood.cpp	
Point_defects.o: Point_defects.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Point_defects.cpp	
Linear_defects.o: Linear_defects.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Linear_defects.cpp	
Strain.o: Strain.cpp
		$(CXX) $(CFLAGS) $(INCLUDES) -c Strain.cpp	


run: main.o IO_data.o Watch.o Hamiltonian.o Psy.o KuboGreenwood.o Point_defects.o Linear_defects.o Strain.o
		$(CXX) main.o IO_data.o Watch.o Hamiltonian.o Psy.o KuboGreenwood.o Point_defects.o Linear_defects.o Strain.o -o run $(LIBS)
		
