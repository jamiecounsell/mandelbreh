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

###Algorithm Details 

####Approach
The algorithm exploits some facts about coprime numbers. We know that:  

- One of the factors of a coprime number must be less than or equal to the square root of that number.
- If one factor is the square root of **n**, we can conclude that the other prime is also the square root of **n**.
- The only factors of a coprime number are itself, 1, and two prime numbers **p** and **q**, which we are looking for.
- All primes greater than 2 must be odd.
- If a number **p** can evenly divide **n** (that is that `n mod p == 0`), then **p** is one of the primes of **n**.
- If we know one of **n**'s factors, **p**, the other factor is `n/p`.

####Algorithm

The algorithm to find **p** and **q**, two prime factors of a coprime **n** is as follows, where:

* **p** is the first factor
* **q** is the second factor
* **n** is the input coprime number
* **np** is the number of processes
* **my_rank** is the rank of the current process

#####Serial Algorithm

1. Set `p = floor(sqrt(n))`
2. If p is even, subtract 1.
3. Check if `p^2 == n`. If so, return `p == q == sqrt(n)`
4. Set `p = p - 2` (jump to next odd number)
5. Check if `n % p == 0`. If so, return `p = p, q = n/p`
6. Jump to step 3

#####Distributed Algorithm

1. Set `p = floor(sqrt(n))`
2. If p is even, subtract 1.
2. For each process, set `p = p - 2*my_rank`
2. Check if `p^2 == n`. If so, return `p == q == sqrt(n)`
3. Set `p = p - 2 * my_rank` (jump `np` odd numbers)
4. Check if `n % p == 0`. If so, return `p = p, q = n/p`
5. Jump to step 3

#####Pseudo Code

```
p = floor(sqrt(n));
p = p - my_rank * 2;

while p is greater than 2 {
    every 1000 loops {
        check_for_message;
        break_if_message;
    }
    if (n % p == 0) {
         found = true;
        notify_others;
        break;
    }
    p = p - 1;
    r++;
}
```

###Results

Using the input co-prime number `220344363094970858831` yielded the following results\*:

|# processes | max_time (s) | speedup | efficiency | %E   |
|------------|----------|---------|------------|------|
| 1          | 915.7197 | 1       | 1          | 100% |
| 2          | 466.7278 | 1.962   | 0.981      | 98%  |
| 4          | 233.5162 | 3.921   | 0.980      | 98%  |
| 8          | 176.2689 | 5.195   | 0.649      | 65%  |
| 16         | 122.8960 | 7.451   | 0.466      | 47%  |
| 32         | 74.80580 | 12.24   | 0.382      | 38%  |
| 64         | 39.70061 | 23.07   | 0.360      | 36%  |

*\*The two prime factors are: (236889251, 4257452468389)*

<img src="http://www.jamiecounsell.me/static/img/mpi/efficiency.png" width=400>

I was unable to find a number that took longer than three hours. I am confident that `71641520761751435455133616475667090434063332228247871795429` would take several hours and would be attainable, but for some reason my processes would always be killed. For the 15 minute result above that took 122.8960s on 16 processes, a breakdown of work is as follows:

```
p   time
--------------
0   122.892453
1   122.896043
2   122.861814
3   122.869866
4   122.894616
5   122.889065
6   122.894171
7   122.859722
8   122.89396
9   122.894442
10  122.885033
11  122.893185
12  122.895691
13  122.87753
14  122.889852
15  122.894295
```
<img src="http://www.jamiecounsell.me/static/img/mpi/cputime.png" width=400>

###Speedup Formula
The answer to problem 2 is attached in hand-written notes.


###Processor Quotes
The following is the answer to problem 3. Our choices are:

* Intel Xeon E5-2630 V3 20MB 8 CORE 2.40GHZ LGA2011 8.00GT/S — $978.98
* Intel Core I7-5820K 3.30GHZ (3.6GHZ Turbo Mode) Six Core 15MB Hyperthreading LGA2011-V3 — $529.99

####Assumptions
- The CPU will be used mostly for floating point computations
- Speed is important
- Money is not a serious issue
- Power consumption is not a serious issue (since money is not)

####Analysis
- Since floating point numbers require little memory relative to other operations, cache size should not be a large concern.
- The E5 is much more power efficient (some sources say 43% more than the I7), but since cost does not matter, this is not important either
- The largest contributor to speed of floating point operations is then the clock rate. The I7 has a much higher clock rate, which would almost surely make it faster.
- Moreover, the I7 comes with an unlocked muliplier, so the CPU can easily be overclocked for greater performance if desired.
- More cores can provide significant speedups, but 8 vs. 6 cores is not as great a difference as 2.4 vs. 3.3 GHz.

```
8 / 6 =      1.333 
3.3 / 2.4 =  1.375
```
- Additionally, the faster clock rate would perform better on programs utilizing a single core, which may occur.

Therefore, I would recommend the Intel Core I7-5820K.