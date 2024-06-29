
#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Naming stuff as variable
code = LinearEquations.cpp
srcdir = src
source = $(srcdir)/matrices.cpp
exe = program.exe
outdir = output
plotdir = plotting
out = $(outdir)/timingoutput.txt

do : $(code) $(source)
	g++ -o $(exe) $<
#	./$(exe) > $(outdir)/output.txt
	./$(exe)
	make plot

output.txt : $(exe) $(source)
	./$< > $@

result.exe : $(code) $(source)
	g++ -o $@ $< 

fresh : 
	mkdir $(outdir) $(plotdir) $(srcdir)



plotscript : 

plot :
	gnuplot gnuplotter
clean :
	rm --force $(exe) $(outdir)/* $(plotdir)/*