You can run it all by just typing "make" into the console

The makefile is clumsily made, but works.
For quickly looking at the results, check the "plotting" folder.

The exercise for the last 1 point was to parallelize the makefile when it runs diagonalization of the N dimension testing.
The "&" makes it run paralelly, and when it is called with -j2 in the "do" recipe it runs using 2 threads.