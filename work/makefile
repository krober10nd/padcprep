FC = mpif90 

libdir2=-L/home2/krober10/DynamicLoadBalance/prep/parmetis-4.0.3/build/Linux-x86_64/libmetis
libdir1=-L/home2/krober10/parmetis/lib/

FFLAGS = -O0  -c -traceback -g -check bounds 
LFLAGS =  -O0  $(libdir1) -lparmetis $(libdir2) -lmetis
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

