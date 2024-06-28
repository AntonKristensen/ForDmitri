
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = roots.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp
plotdir = plotting
out = $(outdir)/roots.txt


do : 
	make clean
	make roots
	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

roots : roots.cpp $(source)
	g++ -o roots.exe $<
#	./roots.exe > $(outdir)/out.txt
	./roots.exe


plot :
	gnuplot gnuplotter

clean :
	rm --force roots.exe $(outdir)/* $(plotdir)/*

