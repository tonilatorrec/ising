FC = gfortran
FLAGS = -no-pie

main.out: main.f random.o
	$(FC) $(FLAGS) main.f random.o -o main.out

random.o: random.f
	$(FC) -c random.f -o random.o

