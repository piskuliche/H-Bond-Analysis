FC=gfortran
FCFLAGS=-Ofast -fopenmp
HOMEPATH=$(PWD)

lb_gmx_inc=/usr2/postdoc/piskulic/Software/gmxfort/libgmxfort-master/bin/include
lb_gmx_lib=/usr2/postdoc/piskulic/Software/gmxfort/libgmxfort-master/bin/lib

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
	chmod 777 bin/*
