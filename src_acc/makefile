VIDEO    = -DVIDEO
CC       = pgc++
CXX=pgc++ 
CPP=pgc++
FLAGS    = -fast -acc -Minfo -ta=tesla,cc35 -O3 $(VIDEO)
CFLAGS   = $(FLAGS)
CXXFLAGS = $(FLAGS)
LDFLAGS  = -acc -ta=tesla,cc35


PROGRAM_NAME=mandelbulb

$(PROGRAM_NAME): main.o print.o timing.o savebmp.o getparams.o 3d.o getcolor.o distance_est.o \
	mandelbulb_dist_est.o raymarching.o renderer.o init3D.o
	#$(CC) $(LDFLAGS) -o $@ $?
	pgc++ -fast -acc -Minfo=accel -ta=tesla,cc35 -O3 -DVIDEO -o mandelbulb *.cc

clean:
	rm -f *.o $(PROGRAM_NAME) $(EXEEXT) *~
