FLAGS=-fopenmp
LAPACK=-llapack
PLOT="set term x11;	
PLOT+=p 'Cpop.out' u 1:4 w l lc 'blue' lw 2 title '|R><R|', 
PLOT+='Cpop.out' u 1:5 w l lc 'red' lw 2 title '|L><L|'"

all: vars.o ranlux.o fmap.o
	gfortran -O3 -o fmap.x $^ $(LAPACK) $(FLAGS)

%.o: %.f95
	gfortran -O3 -c $^ $(FLAGS)

clean:
	@rm -fv *.x
	@rm -fv *.o
	@rm -fv *.out
	@rm -fv *.dat
	@rm -fv *.mod
	@rm -fv *.log

run:
	./fmap.x

plot:
	gnuplot -p -e $(PLOT)

test: all run plot
