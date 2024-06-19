
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = eigenvalues.cpp
srcdir = src
source = $(srcdir)/matrices.cpp $(srcdir)/eigen.cpp
exe = program.exe
outdir = output
plotdir = plotting
out = $(outdir)/eigenoutput.txt

rmax = 10
dr = 0.3

do : 
	make clean
	make eigen
	make hamil
	make -j2 timing
	make plot


#output.txt : $(exe) $(source)
#./$< > $@

#program.exe : $(code) $(source)
#g++ -o $@ $< 

rmaxes := $(shell seq 1 0.1 4.9) $(shell seq 5 1 20)
drs := $(shell seq 1 -0.05 0.1) $(shell seq 0.095 -0.005 0.05)


looptest :
	@for i in $(rmaxes) ; do \
		echo $$i ; \
	done

fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

eigen : eigenvalues.cpp $(source)
	g++ -o eigenvalues.exe $<
	./eigenvalues.exe

hamil : hamiltonning.cpp $(source)
	g++ -o hamiltonning.exe $<
#	This loops over different values of rmax, whose resulting energies are put in output/energies.txt, which is then renamed
	@for i in $(rmaxes) ; do \
		echo rmax is $$i ; \
		echo -n "$$i, " >> output/energies.txt ; \
		./hamiltonning.exe $$i $(dr) ; \
	done
	mv output/energies.txt output/rmaxenergies.txt
#	This loops over different values of dr, whose resulting energies are put in output/energies.txt, which is then renamed
	@for i in $(drs) ; do \
		echo dr is $$i ; \
		echo -n "$$i, " >> output/energies.txt ; \
		./hamiltonning.exe $(rmax) $$i ; \
	done
	mv output/energies.txt output/drenergies.txt
	./hamiltonning.exe  10 0.05

maxdim = 150

timing : diagonaltiming.cpp $(source)
	g++ -o diagonaltiming.exe $<
	@touch output/timing.txt
	$(shell > output/timing.txt)
	@for n in $(shell seq 1 1 $(maxdim)); do \
		./diagonaltiming.exe $$n >> output/timing.txt & \
	done
	@wait
	

plot :
	gnuplot gnuplotter

clean :
	rm --force eigenvalues.exe hamiltonning.exe $(outdir)/* $(plotdir)/*

