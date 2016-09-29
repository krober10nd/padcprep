FC = mpif90 
FFLAGS = -O0 -c -traceback -g -check bounds 
LFLAGS =  -O0
LINK = $(FC)
TARGET =  ../DynExchange.exe 
OBJS =  ../src/MESSENGER.o \
	../prep/PRE.o \
	../prep/PARPREP.o \
        ../src/DRIVER.o \
	
all: $(TARGET)

clean: 
	-rm -rf *.o *.mod ../prep/*.o ../prep/*.mod ../src/*.o ../src/*.mod
$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LFLAGS) 
%.o : $(srcdir)%.F90
	$(FC) $(FFLAGS) $< -o $@

