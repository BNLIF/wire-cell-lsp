MODULES  = givens_rotation householder_transform least_squares lsp nnls qr_decomposition singular_decomposition
CXXFLAGS += -I ../include -Wall -DTEST

all: $(MODULES)

% : %.o
	$(CXX) $(LDFLAGS) $+ -o $@

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~

.PHONY: clean all
