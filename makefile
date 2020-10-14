FLAGS=-fopenmp
LAPACK=-llapack
PLOT="set xrange [-10:2100]; set yrange [-0.1:1.1];
PLOT+=set term x11 size 1000,333; unset key;
PLOT+=p 'exact' u 1:2 w l lc 'black' lw 2, 
PLOT+='exact' u 1:3 w l lc 'black' dt 2 lw 2, 
PLOT+='Cpop.out' u 1:4 w l lc 'blue' lw 2, 
PLOT+='Cpop.out' u 1:5 w l lc 'blue' lw 2 dt 2"

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

plot:
	@gnuplot -p -e $(PLOT)
