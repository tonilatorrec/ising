FC = gfortran
FLAGS = -no-pie -pg 

main.out: main.f 
	$(FC) -o main.out $(FLAGS)