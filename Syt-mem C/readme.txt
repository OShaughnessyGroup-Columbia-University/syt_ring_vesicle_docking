3D Syt+membrane/vesicle simulation by Jie Zhu

To compile main.c under linux: 
	cc -O3 main.c -lm -fopenmp -ffast-math -funroll-all-loops

Running options can be changed by editing the #define... in main.c

To run in parallel, define OPENMP=1

To run with graphics, define OPENGL=1 (requiring GL and glut packages)

To save graphics, define SAVPIC=1 (requiring Devil package)

Parameters are define in paras.ini


