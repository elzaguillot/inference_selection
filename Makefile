TARGET = moran

## the path below must be modified to match the directory of the Eigen library
EIGENPATH = /home/eguillot/Eigen
INCLUDES = -I /usr/include/mysql -I /usr/local/include  -I ${EIGENPATH}
CFLAGS = -O0 -g -pg -std=c++11 ${INCLUDES}

# library flags go here
LFLAGS = -O0 -pg -g

CC = g++

SRCDIR = src

# specifying source files
SRC = src/main.cpp src/moran.cpp src/moranMatrix.cpp src/wfisher.cpp src/utils.cpp
SRC2 = src/compare_main.cpp src/moran.cpp src/moranMatrix.cpp src/wfisher.cpp src/utils.cpp src/wrightFisher.cpp 
SRC3 = src/vscompute_main.cpp src/moran.cpp src/moranMatrix.cpp src/wfisher.cpp src/utils.cpp
SRC4 = src/rmain.cpp src/moran.cpp src/moranMatrix.cpp src/wfisher.cpp src/utils.cpp

# suffix .cpp in $SRC will be substitued to .o
OBJ = $(SRC:.cpp=.o)
OBJ2 = $(SRC2:.cpp=.o)


# rule for creating the final binary. $OBJ is a list of *.o files obtained by compiling *.cpp files

$(TARGET): $(OBJ)
	@echo Linking $@
	@$(CC) $(LFLAGS) $(OBJ) -o $@
	@echo Build Complete

# rule to compile *.o object files from *.cpp files
.cpp.o: $<
	@echo Compiling $<
	@$(CC) -c $(CFLAGS) $< -o $@

compare: $(OBJ2)
	@echo Linking $@
	@$(CC) $(LFLAGS) $(OBJ2) -o $@
	@echo Build Complete

vscompute: $(OBJ3)
	@echo Linking $@
	@$(CC) $(LFLAGS) $(OBJ3) -o $@
	@echo Build Complete

rlib: $(SRC4)
	cd src/
	export PKG_LIBS='`Rscript -e "Rcpp:::LdFlags()"` -O2 -std=c++11 -I /home/eguillot/Eigen'
	export PKG_CXXFLAGS='`Rscript -e "Rcpp:::CxxFlags()"` -O2 -std=c++11 -I /home/eguillot/Eigen'
	@R CMD SHLIB $(SRC4) #moranMatrix.cpp rmain.cpp wfisher.cpp utils.cpp



# rule for 'make clean' command
clean:
	@rm -f $(OBJ) $(TARGET) $(OBJ2) $(TARGET2)
	@echo All object files and binaries removed
