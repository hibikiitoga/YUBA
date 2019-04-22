CXX=clang++-mp-7.0
SRCS=main.cpp
CXXFLAGS=-std=c++2a -stdlib=libc++ -I/opt/local/include/ -pthread -I./eigen/
EXE_FILE=YUBA
OBJS=$(SRCS:.cpp=.o)
DEPS=$(SRCS:.cpp=.d)
-include $(DEPS)

.PHONY: debug
debug: CXXFLAGS+=-O2 -Wall -Wextra -Wshadow -pedantic
debug: build

.PHONY: release
release: CXXFLAGS+=-O2 -DIn_RELEASE=true
release: build

.PHONY: build
build: $(OBJS)
	$(CXX) -o $(EXE_FILE) $(OBJS) Chaperone/heap.o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(SETTINGS) $<

.PHONY: clean
clean:
	rm -f $(EXE_FILE)* $(OBJS)
