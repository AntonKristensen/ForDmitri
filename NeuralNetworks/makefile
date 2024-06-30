
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = neural.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp
plotdir = plotting
out = $(outdir)/neural.txt


do : 
	make clean
	make neural
#	make plot



fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

neural : neural.cpp $(source)
	g++ -o neural.exe $<
#	./neural.exe > $(outdir)/out.txt
	./neural.exe


plot :
	gnuplot gnuplotter

clean :
	rm --force neural.exe $(outdir)/* $(plotdir)/*

