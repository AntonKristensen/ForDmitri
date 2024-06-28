
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = minimization.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp
plotdir = plotting
out = $(outdir)/minimization.txt


do : 
	make clean
	make minimization
	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

minimization : minimization.cpp $(source)
	g++ -o minimization.exe $<
	./minimization.exe > $(outdir)/out.txt
#	./minimization.exe


plot :
	gnuplot gnuplotter

clean :
	rm --force minimization.exe $(outdir)/* $(plotdir)/*

