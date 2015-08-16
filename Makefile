CXXFLAGS=-g	
INCDIR=/Users/nishita/Code/include
LIBDIR=/Users/nishita/Code/lib
LDFLAGS1 := $(shell root-config --ldflags --glibs) 

atlasmet: atlasmet.cc
	$(CC) $(CXXFLAGS) -I $(INCDIR) atlasmet.cc -o atlasmet.x \
	-lstdc++ -L /Users/nishita/Code/Py8-trunk/lib/archive -lpythia8 -llhapdfdummy

stau: 	stau.cc
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -o $@.exe $@.cc \
	-L$(LIBDIR) -lpythia8 -llhapdfdummy $(LDFLAGS1)
