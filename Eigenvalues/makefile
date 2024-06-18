
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Naming stuff as variable
code = eigenvalues.cpp
srcdir = src
source = $(srcdir)/matrices.cpp $(srcdir)/eigen.cpp
exe = program.exe
outdir = output
plotdir = plotting
out = $(outdir)/eigenoutput.txt


output.txt : $(exe) $(source)
	./$< > $@

program.exe : $(code) $(source)
	g++ -o $@ $< 

fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

do : $(code) $(source)
	g++ -o $(exe) $<
	./$(exe)
	make plot

plotscript : 

plot :
	gnuplot gnuplotter
clean :
	rm --force $(exe) $(outdir)/* $(plotdir)/*