
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = spines.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp
plotdir = plotting
out = $(outdir)/spines.txt


do : 
	make clean
	make integrals
	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

integrals : integrals.cpp $(source)
	g++ -o integrals.exe $<
#	./integrals.exe
	./integrals.exe > $(outdir)/out.txt


plot :
	gnuplot gnuplotter

clean :
	rm --force integrals.exe $(outdir)/* $(plotdir)/*

