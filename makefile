FLAGS=-fopenmp
LAPACK=-llapack
PLOT="set xrange [-100:2300];
PLOT+=set term x11 size 1500,500;	
PLOT+=p 'exact' u 1:2 w l lc 'black' lw 2 title 'Exact', 
PLOT+='exact' u 1:3 w l lc 'black' dt 2 lw 2 notitle, 
PLOT+='Cpop.out' u 1:4 w l lc 'blue' lw 2 dt 2 notitle, 
PLOT+='Cpop.out' u 1:5 w l lc 'blue' lw 2 title 'LSCI/II', 
PLOT+='Cimp.out' u 1:4 w l lc 'red' lw 2 dt 2 notitle, 
PLOT+='Cimp.out' u 1:5 w l lc 'red' lw 2 title 'Improved'"

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
	gnuplot -p -e $(PLOT)
