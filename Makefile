FC=gfortran
FCFLAGS=-Ofast -fopenmp -std=f2008
HOMEPATH=$(PWD)

lb_gmx_inc=/home/piskuliche/Software/libgmxfort/bin/include
lb_gmx_lib=/home/piskuliche/Software/libgmxfort/bin/lib

all: hba hyd

hba: src/fortran/funcs.f90 src/fortran/hbond-finder.f90 
	mkdir -p bin/
	touch bin/test
	rm bin/*
	test -f module/path.include || echo "module use $(HOMEPATH)/module" >> ~/.bash_profile
	test -f module/path.include || source ~/.bash_profile
	@echo "prepend_path('PATH', '$(HOMEPATH)/bin')" > module/path.include
	cat module/hba_header module/path.include > module/hba.lua
	$(FC) $(FCFLAGS) -I $(lb_gmx_inc) -L $(lb_gmx_lib) -lgmxfort -o bin/hba src/fortran/funcs.f90  src/fortran/hbond-finder.f90
	ln -s $(HOMEPATH)/src/python/mark_acceptors.py bin/
	ln -s $(HOMEPATH)/src/python/pull_atoms.py bin/
	ln -s $(HOMEPATH)/src/python/generate_atom_map.py bin/
	ln -s $(HOMEPATH)/src/python/setup_analysis.py bin/
	chmod 777 bin/*

hyd: src/fortran/hydration-shell.f90
	$(FC) $(FCFLAGS) -I $(lb_gmx_inc) -L $(lb_gmx_lib) -lgmxfort -o bin/hyd src/fortran/funcs.f90 src/fortran/hydr_module.f90 src/fortran/hydration-shell.f90
	chmod 777 bin/*
