CC = g++
OBJ = integral.o iso7.o
%.o: %.c
	$(CC) -c $<
nuc.exe: $(OBJ)
		$(CC) -o nuc $(OBJ)
clean:
	rm -f *.o nuc.exe