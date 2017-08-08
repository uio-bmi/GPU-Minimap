# Very simple makefile for now

CC=gcc

CFLAGS=  -O5


minimap: minimap.c minimap.h kseq.h hash.c sort.c LIS.c
	$(CC) $(CFLAGS) -o minimap minimap.c hash.c sort.c LIS.c

debug: minimap.c minimap.h kseq.h hash.c sort.c LIS.c
	$(CC) $(CFLAGS) -g -o minimap minimap.c hash.c sort.c LIS.c

profile: minimap.c minimap.h kseq.h hash.c sort.c LIS.c
	$(CC) $(CFLAGS) -g -pg -o minimap_profile minimap.c hash.c sort.c LIS.c

clean:
	$ rm minimap minimap_profile *.gcno *.gch *.gcov
