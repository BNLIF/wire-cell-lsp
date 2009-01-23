CXXFLAGS += -I ./include

all: givens_rotation.test householder_transform.test singular_decomposition.test least_squares.test qr_decomposition.test lsp nnls.test

givens_rotation.test: givens_rotation.cpp ./include/lsp/givens_rotation.h
	$(CXX) $(CXXFLAGS) -o givens_rotation.test -DTEST givens_rotation.cpp

householder_transform.test: householder_transform.cpp ./include/lsp/householder_transform.h ./include/lsp/singular_decomposition.h
	$(CXX) $(CXXFLAGS)  -o householder_transform.test -DTEST householder_transform.cpp

singular_decomposition.test: singular_decomposition.cpp ./include/lsp/singular_decomposition.h ./include/lsp/householder_transform.h ./include/lsp/qr_decomposition.h
	$(CXX) $(CXXFLAGS)  -o singular_decomposition.test -DTEST singular_decomposition.cpp

qr_decomposition.test: qr_decomposition.cpp ./include/lsp/qr_decomposition.h ./include/lsp/givens_rotation.h
	$(CXX) $(CXXFLAGS)  -o qr_decomposition.test -DTEST qr_decomposition.cpp

least_squares.test: least_squares.cpp ./include/lsp/least_squares.h ./include/lsp/singular_decomposition.h
	$(CXX) $(CXXFLAGS)  -o least_squares.test -DTEST least_squares.cpp

nnls.test: nnls.cpp ./include/lsp/nnls.h ./include/lsp/least_squares.h
	$(CXX) $(CXXFLAGS)  -o nnls.test -DTEST nnls.cpp

lsp: lsp.cpp
	$(CXX) $(CXXFLAGS)  -o lsp -DNDEBUG lsp.cpp

clean:
	rm -f *.o
	rm -f *~
	rm -f *.test

.PHONY: clean all
