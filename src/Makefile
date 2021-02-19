OBJS = domain_t.o example_set.o fuzzy_t.o genetcode.o multipoblacion.o poblacionbinaria.o poblacionentera.o poblacionreal.o Relaciones.o ruleset.o variable_t.o vectordouble.o VectorVar.o prueba.o operaciones.o patrones.o


all: $(OBJS)
	g++ -std=c++14 -fopenmp -fpermissive -O2 -o WyM_Ver $(OBJS)

%.o : %.cpp
	g++ -std=c++14 -fopenmp -fpermissive -O2 -c $< -o $@

domain_t.o: domain_t.cpp
example_set.o: example_set.cpp
fuzzy_t.o: fuzzy_t.cpp
genetcode.o: genetcode.cpp
multipoblacion.o: multipoblacion.cpp
poblacionbinaria.o: poblacionbinaria.cpp
poblacionentera.o: poblacionentera.cpp
poblacionreal.o: poblacionreal.cpp
Relaciones.o: Relaciones.cpp
ruleset.o: ruleset.cpp
variable_t.o: variable_t.cpp
vectordouble.o: vectordouble.cpp
VectorVar.o: VectorVar.cpp
prueba.o: prueba.cpp
operaciones.o: operaciones.cpp
patrones.o: patrones.cpp

clean:
	 rm -f $(OBJS) WyM_Ver
