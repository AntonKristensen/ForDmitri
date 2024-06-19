
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = least.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp $(srcdir)/eigen.cpp
plotdir = plotting
out = $(outdir)/least.txt


do : 
	make clean
	make least
	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

least : least.cpp $(source)
	g++ -o least.exe $<
	./least.exe


plot :
	gnuplot gnuplotter

clean :
	rm --force eigenvalues.exe hamiltonning.exe $(outdir)/* $(plotdir)/*

