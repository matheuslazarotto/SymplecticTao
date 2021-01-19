# Folder with '.hpp' files
INCLUDE_DIR = ./include

# Folder with '.cpp' files
VPATH = ./src

# Linked libraries
LINKFLAGS = -lm

# Compilation flags
# -w : supresses all warnings
# -Wall : turn on compilar warning flags
CXXFLAGS = -g -Wall -I$(INCLUDE_DIR)

# Dependents  files
DEPENDENCIES = main.cpp ode_sys_solv.cpp auxf.cpp

# Output (executable) 
TARGETS = symplectic_tao

build: $(TARGETS)

$(TARGETS) : $(DEPENDENCIES)
	$(CXX) $(CXXFLAGS) -o $(TARGETS) $^ $(LINKFLAGS)

.PHONY: clean
clean:
	-rm -f $(TARGETS)