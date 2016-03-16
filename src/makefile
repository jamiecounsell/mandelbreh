#CFLAGS=-g -O2 -Wall -pg -DMONITOR

all: ser acc

ser: VIDEO    = -DVIDEO
ser: CXX      =  g++  
ser: FLAGS    = -O3 $(VIDEO) -Wall
ser: CFLAGS   = $(FLAGS)
ser: CXXFLAGS = $(FLAGS)
ser: LDFLAGS  = -lm

acc: VIDEO    = -DVIDEO
acc: CC       = pgc++
acc: CXX=pgc++ 
acc: CPP=pgc++
acc: FLAGS    = -fast -acc -Minfo=accel -Minline -ta=tesla,cc35 -O3 $(VIDEO)
acc: CFLAGS   = $(FLAGS)
acc: CXXFLAGS = $(FLAGS)
acc: LDFLAGS  = -acc -Minfo=accel -Minline -ta=tesla,cc35

ser: PROGRAM_NAME=mandelbulb_ser
acc: PROGRAM_NAME=mandelbulb_acc

OBJS = 
MP_OBJS = main_mp.o readwriteppm.o compimage_mp.o
ACC_OBJS = main_acc.o readwriteppm.o compimage_acc.o


#$(PROGRAM_NAME): $(OBJS)
ser: $(OBJS)
	$(CC) $(LDFLAGS) -o mandelbulb_ser $?

acc: main_acc.o readwriteppm.o compimage_acc.o
	$(CC) $(LDFLAGS) -o genimg_acc $?

clean:
  rm -f *.o mandelbulb_acc mandelbulb_ser $(EXEEXT) *~