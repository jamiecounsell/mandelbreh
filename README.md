##4F03 Project

###Authors  
  
| Name           | Student Number | Email                | Website                                     |
|:---------------|:---------------|:---------------------|:--------------------------------------------|
| Jamie Counsell | 1054209        | counsej@mcmaster.ca | [jamiecounsell.me](http://jamiecounsell.me/) |
| James Priebe | 1135001        | priebejp@mcmaster.ca |  |

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

* **mandelbulb** - Generate a Mandelbulb fractal
* **mandelbox** - Generate a Mandelbox fractal
* **boxserial** - Generate a Mandelbox fractal using a serial implementation
* **bulbserial** - Generate a Mandelbulb fractal using a serial implementation

Then run the command with the optional runtime flags:

```
$ ./mandel[box, bulb, bulb_serial, box_serial] [-f n] [-v] [-n]
```
where:

* **f** - Instruct the program to generate `n` frames. Default is 1 frame.  
* **v** - Instruct the program to generate a video when it is complete (calling `genvideo.sh`)
* **n** - Is a less verbose option. Use this to silence more output.
