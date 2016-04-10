##4F03 Project

###Authors  
  
| Name           | Student Number | Email                | Website                                      |
|:---------------|:---------------|:---------------------|:---------------------------------------------|
| Jamie Counsell | 1054209        | counsej@mcmaster.ca  | [jamiecounsell.me](http://jamiecounsell.me/) |
| James Priebe   | 1135001        | priebejp@mcmaster.ca |                                              |

###Description
This C program renders a mandelbulb or mandelbox 3D fractal with optimizations in OpenACC.

###Dependencies
- `pgcc` or compiler with support for OpenACC

###Installation
- Install the dependencies
- Clone the repository

###Operation
To run the program, first run make:

```
$ make clean
$ make [type]
```
Where `type` is one of:

* **mandelbulb** - Compute a Mandelbulb fractal
* **mandelbox** - Compute a Mandelbox fractal
* **boxserial** - Compute a Mandelbox fractal using a serial implementation
* **bulbserial** - Compute a Mandelbulb fractal using a serial implementation

Then run the command with the optional runtime flags:

```
$ ./mandel[box, bulb, bulb_serial, box_serial] params.dat [-f n] [-v]
```
where:

* **params.dat** is a file containing the Mandelbulb or Mandelbox parameters. 
* **f** - Instruct the program to generate `n` frames. Default is 1 frame.  
* **v** - Instruct the program to generate a video when it is complete (calling `genvideo.sh`)

For example, to generate a 7200 frame video at 30FPS (4 minutes) uwisng the mandelbulb:

```
$ ./mandelbulb bulb_params.dat -f 7200 -v
```

###Speedups

For the first frame of the submitted video, the following times were recorded:

|Server|OpenACC|time      |
|------|-------|----------|
|tesla |NO     |108.16456s|
|tesla |YES    |1.236812s |

###Parallelization
The only region that was parallelized was the nested loop in `renderer.cc`. This loop is the program's largest bottleneck and also supports parallelization quite intuitively. OpenACC pragmas were used to identify the region as an OpenACC compute region, as well as transfer the data to the device from the host. The outer loop was explicitly marked as parallel, and other optimizations were left up to PGCC.  

Since frame parameters were not generated asynchronously, no parallelization was done to compute more than one frame at a time.

###Frame Generation
Frames are generated sequentially from an array of `CameraParams` structures. The first image generated is always the same as what is identified in the input parameters. This ensures that the assignment requirements can be properly met with the given `paramsBulb.dat` file. After the first image, the camera rotates around the fractal, slowly decreasing its position in the `z` axis from `1` to `-1` across 7200 frames. The position is computed as follows:

* the `x` coordinate is `cos(frame_number/500)`
* the `y` coordinate is `sin(frame_number/500)`
* the `z` coordinate is `1 - (frame_number/3600)`

This will guide the camera around the object in a circular motion along the `x, y` plane such that it will complete one full rotation every `500*pi` frames. The z value decreases individually from `1` to `-1` between frames `0` and `7200`, respectively. This creates a sort of _spring_ path, showcasing all sides of the fractal.

Each iteration, `init3D` is called again to ensure the camera is still facing the center point at `(0,0,0)`.


###Final Result


###Source Code

See attached.
