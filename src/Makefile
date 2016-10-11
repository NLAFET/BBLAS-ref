include ../make.inc

DIR_DEPS = ../include
INC	 = -I $(DIR_DEPS)

#Source code files
BBLAS_SRC= zgemm_batch.c xerbla_batch.c

#Object files
BBLAS_OBJ=$(BBLAS_SRC:.c=.o)

#Executable
EXE=test

# Targets to be created
all: $(BBLAS_OBJ)


#Object files creation
.c.o:
	$(CC) $(CFLAGS) -fopenmp $(INC) $<  -o $@

clean:
	rm *.o *~
	rm $(EXE)
