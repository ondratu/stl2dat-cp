SRC = $(wildcard *.cxx)
OBJ = $(SRC:%.cxx=%.o)

all: stl2dat

%.o: %.cxx
	$(CXX) -MMD  $(CFLAGS) -c $< -o $@

stl2dat: $(OBJ) $(wildcard *.h)
	$(CXX) -o stl2dat $(OBJ)
