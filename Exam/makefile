# A decent amount of stuff in here is just relics from previous iterations. Sorry that I haven't cleaned it up :(


#dish : ingredients
#⊣recipe

# $@ is the file name of the target of the rule.
# $< is the name of the first prerequisite
# $^ are the names of all the prerequisites

# Passing a variable from commandline into the make program
# "make foo=bar target"
#Then in your makefile you can refer to $(foo)

# Naming stuff as variable
code = exam.cpp
srcdir = src
outdir = output
source = $(srcdir)/matrices.cpp 
plotdir = plotting
out = $(outdir)/exam.txt


do : 
	make clean
	make -j6 scaling
	make plot

ns := $(shell seq 10 10 500) $(shell seq 550 50 1000) $(shell seq 1250 250 3000)
scaling : exam.cpp
	g++ -o exam.exe $<
	@for i in $(ns) ; do \
		echo "N is $$i/3000, "  ; \
		./exam.exe $$i >> $(outdir)/out.txt; \
	done


fresh : 
	mkdir -p $(outdir) $(plotdir) $(srcdir)

exam : exam.cpp $(source)
	g++ -o exam.exe $<
#	./exam.exe > $(outdir)/out.txt
	./exam.exe 100


plot :
	gnuplot gnuplotter

clean :
	rm --force exam.exe $(outdir)/* $(plotdir)/*

