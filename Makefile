FC = gfortran
FLAGS = -no-pie

main.out: main.f 
	$(FC) -o main.out $(FLAGS) main.f