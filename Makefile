CXXFLAGS += -I ./include

all: givens_rotation.test householder_transform.test singular_decomposition.test least_squares.test qr_decomposition.test lsp nnls.test

givens_rotation.test: givens_rotation.cpp givens_rotation.h
	$(CXX) $(CXXFLAGS) -o givens_rotation.test -DTEST givens_rotation.cpp

householder_transform.test: householder_transform.cpp householder_transform.h singular_decomposition.h
	$(CXX) $(CXXFLAGS)  -o householder_transform.test -DTEST householder_transform.cpp

singular_decomposition.test: singular_decomposition.cpp singular_decomposition.h householder_transform.h qr_decomposition.h
	$(CXX) $(CXXFLAGS)  -o singular_decomposition.test -DTEST singular_decomposition.cpp

qr_decomposition.test: qr_decomposition.cpp qr_decomposition.h givens_rotation.h
	$(CXX) $(CXXFLAGS)  -o qr_decomposition.test -DTEST qr_decomposition.cpp

least_squares.test: least_squares.cpp least_squares.h singular_decomposition.h
	$(CXX) $(CXXFLAGS)  -o least_squares.test -DTEST least_squares.cpp

nnls.test: nnls.cpp nnls.h least_squares.h
	$(CXX) $(CXXFLAGS)  -o nnls.test -DTEST nnls.cpp

lsp: lsp.cpp
	$(CXX) $(CXXFLAGS)  -o lsp -DNDEBUG lsp.cpp

clean:
	rm -f *.o
	rm -f *~
	rm -f *.test

.PHONY: clean all
