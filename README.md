# pace-2022

***

**DiVerSeS** - a solver for the directed feedback vertex set problem, written as an entry to the [PACE 2022 challenge](https://pacechallenge.org/).
This repository contains code of an exact solver and a heuristic solver.

To use exact solver, please use '-\-track=exact' option. In order to use heuristic solver, please use '-\-track=heuristic' (if not specified, then by default a heuristic solver will be run).

The heuristic solver was optimized to provide best results within 10 minutes on optil.io platform. 
If you would like to use the solver in other 'settings', do not hesitate to let me know, I might be able to help you to choose best possible configuration to meet your needs and make DiVerSeS work better in your specific context. Other 'settings' mean e.g. running solver for time limit different than the default 590 seconds, running the solver for some specific graph class (e.g. large sparse random graphs or real-world networks), trying first to produce kernel of smallest possible size using advanced data reduction rules - which may be time consuming for very large graphs and is not used by default - etc. If you have any questions or suggestions, just let me know.

***

**Requirements**:

CMake VERSION 3.10.2 or higher<br>
c++ 17 or higher<br>

For best performance of DiVerSeS you will also need to have executable binaries of the following solvers:<br>
An executable binary of a vertex cover solver WeGotYouCovered (solver available [here](https://github.com/KarlsruheMIS/pace-2019), code used - unchanged - can be found in directory WGYC), named "vc-solver" should be placed in the same directory as DiVerSeS solver.<br>
An executable binary of a hitting set solver FINDMINHS (solver available [here](https://github.com/felerius/findminhs), code used - slightly modified - can be found in directory FINDMINHS), named "findminhs" should be placed in the same directory as DiVerSeS solver.<br>

You may use DiVerSeS without those solvers, but it may be slower on the exact track (when run without WeGotYouCovered solver it will be roughly 2-3 times slower on average when run on public instances of PACE 2022 challenge, when run without FINDMINHS solver it will be loads slower on 'not-kernelizable' graphs with 'low' percentage of bidirectional arcs). The performance of a heuristic solver DiVerSeS should not be visibly affected by not using the two solvers. (please also note that DiVerSeS was not highly tested for use without the two solvers). <br>

***

**Installation**:


Use cmake to obtain a binary file, e.g. in linux in the main directory (pace-2022) you can use the following commands:

mkdir build<br>
cd build<br>
cmake ..<br>
make <br>
<br>
After this, the executable file named "DiVerSeS" should be in the "build" directory.

Optional: build solvers mentioned in 'Requirements' section and place properly named binaries in the directory where DiVerSeS can be found (in the example this is the "build" directory).

***

**Usage:**

Given a graph in a file example_input.txt, you can run DiVerSeS e.g. in the following way
 
./DiVerSeS < example_input.txt > example_output.out 2>example_logs.err

DiVerSeS (version from heuristic track, on main branch) will run for 590 seconds (default) or until it receives SIGTERM signal (e.g. using "kill -s SIGTERM $pid" command in linux).<br>
The exact version of DiVerSeS will terminate for itself and does not need SIGTERM to be sent to the process. Please do not send the signal to the process if you use an exact algorithm.
 
You can specify the time limit (in milliseconds) for the heuristic solver DiVerSeS with '-\-time-limit' option. Please bear in mind that this limit may not be observed if you run the solver on large graphs with relatively small time limit. Changing just time limit to value other than the default (and not adapting properly the workflow of DiVerSeS) may not be an optimal way to produce best results within the specified time limit (especailly for large graphs) - the heuristic solver was optimized to provide best results within 10 minutes on optil.io platform.

To use exact solver, please specify '-\-track=exact' option. You can disable logs using '-\-quiet=true' option. You can also provide a path to the file with input data using '-\-file=some_path' option (if not provided, then solver will read data from standard input).

<br>

Exemplary usage with additional parameters:

./DiVerSeS -\-time-limit=60000 -\-file=example_input.txt > example_output_heuristic.out

./DiVerSeS -\-quiet=true -\-track=exact < example_input.txt > example_output_exact.out
