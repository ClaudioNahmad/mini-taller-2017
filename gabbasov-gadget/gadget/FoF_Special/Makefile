
#OPT1   =  -DVELOCITYDISPERSION
#OPT2   =  -DPERIODIC

EXEC   = FoF_Special

OBJS   = main.o io_input.o allvars.o  allocate.o peano.o ngbtree.o \
	 nrsrc/nrutil.o nrsrc/sort2_int.o nrsrc/sort_int.o nrsrc/sort2_fltint.o nrsrc/iindexx.o


INCL   = allvars.h   nrsrc/nrutil.h proto.h Makefile

OPTIONS =  $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)

CFLAGS =  -O2 -Wall  $(OPTIONS)

#SYSTYPE="Regatta"

ifeq ($(SYSTYPE),"Regatta")
CFLAGS =  -O2 -q64  $(OPTIONS)
endif



LIBS   =  -lm 

CC     =  cc


$(EXEC): $(OBJS) 
	cc $(CFLAGS) $(OBJS)  -o $(EXEC)  $(LIBS)

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
