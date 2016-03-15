
CC       =  g++
FLAGS    = -O3 -Wall
CFLAGS   = $(FLAGS)
CXXFLAGS = $(FLAGS)
LDFLAGS  = -lm 

PROGRAM_NAME=mandelbulb

$(PROGRAM_NAME): main.o print.o timing.o savebmp.o getparams.o 3d.o getcolor.o distance_est.o \
	mandelbulb_dist_est.o raymarching.o renderer.o init3D.o
	$(CC) -o $@ $? $(CFLAGS) $(LDFLAGS)


clean:
	rm *.o $(PROGRAM_NAME) $(EXEEXT) *~
