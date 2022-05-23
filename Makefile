FC = gfortran
FLAGS = -no-pie -pg 

main: main.f 
	$(FC) -o main $(FLAGS)