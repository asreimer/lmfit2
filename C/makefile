######### PHONY ##########

.PHONY: all clean testing 

######### VARIABLES #########

CC = gcc

#CFLAGS = -g -pedantic -Wall -D_GNU_SOURCE $(INCLUDES)
 
CFLAGS = -O3 -pedantic -Wall -D_GNU_SOURCE $(INCLUDES) 

INCLUDES = -Iinclude/ -I$(IPATH)/superdarn/ -I$(IPATH)/base/ -I$(IPATH)/general/

OBJECTDIR = objects/

SRCDIR = src/

LIBDIR = lib/

BINDIR = bin/

STATIC_LFLAGS = -Wl,-Bstatic $(LIBDIR)liblmfit2.a

DYNAMIC_LFLAGS = -Wl,-Bdynamic -lpthread -lm

MFLAGS = -L$(LIBDIR) 

LIBS =  $(MFLAGS) $(STATIC_LFLAGS) $(DYNAMIC_LFLAGS)

######### WHAT TO COMPILE #########

all: lmfit2 make_lmfit2


######### LMFIT2 LIBRARY ##########
LMFIT2_TARGET = liblmfit2.a

LMFIT2_SRC = $(shell echo src/*.c)

LMFIT2_OBJ = $(LMFIT2_SRC:.c=.o)

lmfit2:$(LMFIT2_OBJ)
	mkdir -p $(LIBDIR)
	ar -rc $(LMFIT2_TARGET) $(LMFIT2_OBJ)
	ranlib $(LMFIT2_TARGET)
	mv $(LMFIT2_TARGET) $(LIBDIR)


######### TESTING ##########
TESTING = testing_lmfit2

TESTING_SRC = $(shell echo testing/*.c)

TESTING_OBJ = $(TESTING_SRC:.c=.o)

testing: $(TESTING_OBJ)
		$(CC) $(CFLAGS) $(INLUDES) -o $(TESTING) $(TESTING_OBJ) $(LIBS)

######### MAKE_FIT #########
MAKEFIT = make_lmfit2

MAKEFIT_SRC = $(shell echo make_lmfit2.0/*.c)

MAKEFIT_OBJ = $(MAKEFIT_SRC:.c=.o)

MK_STATIC_LFLAGS = -Wl,-Bstatic $(LIBDIR)liblmfit2.a -loldraw.1 -loldfit.1 -lcfit.1 -lrscan.1 -lradar.1 -lraw.1 -lfit.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 

MK_DYNAMIC_LFLAGS = -Wl,-Bdynamic -lpthread -lm -lz

MFLAGS = -L$(LIBDIR) -L$(LIBPATH)

MK_LIBS =  $(MFLAGS) $(MK_STATIC_LFLAGS) $(MK_DYNAMIC_LFLAGS)

make_lmfit2: $(MAKEFIT_OBJ)
		     $(CC) $(CFLAGS) $(INLUDES) -o $(MAKEFIT) $(MAKEFIT_OBJ) $(MK_LIBS)
		     mkdir -p $(BINDIR)
		     mv $(MAKEFIT) $(BINDIR)


######### CLEAN #########
clean: 
	rm -rf *.h~ *.c~ *.o *.o~ *.a~ *.exe *.exe.stackdump Makefile~ fullrangeinfo.txt phases.txt
	rm -rf *.*.gch *.txt~  $(TESTING) $(MAKEFIT) $(LIBDIR)$(LMFIT2_TARGET) $(LMFIT2_OBJ) $(TESTING_OBJ) $(MAKEFIT_OBJ)




