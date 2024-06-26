
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = montecarlo.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp
plotdir = plotting
out = $(outdir)/montecarlo.txt


do : 
	make clean
	make montecarlo
	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

montecarlo : montecarlo.cpp $(source)
	g++ -o montecarlo.exe $<
	./montecarlo.exe > $(outdir)/out.txt
#	./montecarlo.exe


plot :
	gnuplot gnuplotter

clean :
	rm --force montecarlo.exe $(outdir)/* $(plotdir)/*

