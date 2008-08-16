all: givens_rotation.test householder_transform.test singular_decomposition.test least_squares.test

givens_rotation.test: givens_rotation.cpp givens_rotation.h
	$(CXX) -o givens_rotation.test -DTEST givens_rotation.cpp

householder_transform.test: householder_transform.cpp householder_transform.h
	$(CXX) -o householder_transform.test -DTEST householder_transform.cpp

singular_decomposition.test: singular_decomposition.cpp singular_decomposition.h
	$(CXX) -o singular_decomposition.test -DTEST singular_decomposition.cpp

least_squares.test: least_squares.cpp least_squares.h
	$(CXX) -o least_squares.test -DTEST least_squares.cpp

clean:
	rm -f *.o
	rm -f *~

.PHONY: clean all
