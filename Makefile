CC = g++
OBJ = BD.o iso7.o screen5.o
%.o: %.c
	$(CC) -c $<
nuc.exe: $(OBJ)
		$(CC) -o nuc.exe $(OBJ)
clean:
	rm -f *.o nuc.exe
