all: permute.so merge_particle_stilt_files

permute.so: permute_c.c
	gcc -O3 -I/usr/share/R/include -rdynamic -fPIC -shared permute_c.c -o permute.so

merge_particle_stilt_files: merge_particle_stilt_files.cpp
	g++ -O3 -Wall merge_particle_stilt_files.cpp -o merge_particle_stilt_files
