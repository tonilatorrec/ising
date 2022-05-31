FC = gfortran
FLAGS = -no-pie

main.out: main.f random.o thermo.o
	$(FC) $(FLAGS) main.f random.o -o main.out

random.o: random.f
	$(FC) -c random.f -o random.o

thermo.o: thermo.f
	$(FC) -c thermo.f -o thermo.o