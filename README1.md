<div align="center">
  	<h1>
    	Artifact Evaluation
  	</h1>
</div>

#  Artifact Requirements

##  Resource requirements

 All  experiments herein are run on Ubuntu 20.04.6 LTS in virtual machine 
with CPU 12th Gen Intel Core i9-12900K × 16  and RAM 15.6 GB. The compared tool CORA is run in Matlab R2020a. Here are the **third-party libraries** required for BdryReach.
| Library | Website | Version |
| --- | --- | --- |
| git | [https://git-scm.com/](https://git-scm.com/) | 2.25.1 |
| Cmake | [https://cmake.org/](https://cmake.org/) | latest |
| Eigen | [http://eigen.tuxfamily.org/index.php?title=Main Page](http://eigen.tuxfamily.org/index.php?title=Main%20Page) | 3.34 |
| Python | [https://www.python.org/](https://www.python.org/) | 3.8 |
| Capd | [http://capd.ii.uj.edu.pl/](http://capd.ii.uj.edu.pl/) | latest |
| boost | [https://www.boost.org/](https://www.boost.org/) | 1.67.0.0 |
| GLPK | [https://www.gnu.org/software/glpk/](https://www.gnu.org/software/glpk/) | 4.65-2 |

##  Evaluation runtime

We provided 8 automation scripts to run all the examples in our paper for the artifact evaluation. The examples files are located in **./examples**, which contains all the examples showcased in Table 3-5. It has tow folders, **./examples/BdryReach** and **./examples/CORA**, which include examples run by BdryReach (our tool) and CORA (compared tool) respectively.

In the directory **./examples/BdryReach**, it has 4 automation scripts corresponding to examples in Table 3, Table 4, examples with short time horizon in Table 5 and examples with long time horizon in Table 5.
Here are the runtime of these 4 scripts, which are tested in the environment stated in 1.1.
| Automation scripts | Runtime (seconds) |
| --- | --- |
| runTable3_Examples.sh | 425 | 
| runTable4_Examples.sh | 1453 | 
|runTable5_Examples_shortTime.sh | 527 | 
| runTable5_Examples_longtTime.sh  | 840 | 

In the directory **./examples/CORA**, it has 4 automation scripts corresponding to examples in Table 3, Table 4, examples with short time horizon in Table 5 and examples with long time horizon in Table 5. 
Here are the runtime of these 4 scripts, which are tested in the environment stated in 1.1.

| Automation scripts | Runtime (seconds) |
| --- | --- |
| runTable3_Examples.m | 2800 | 
| runTable4_Examples.m | 409 | 
|runTable5_Examples_shortTime.m | 4484 | 
| runTable5_Examples_longtTime.m  | 390 | 
# Structure and Content

- [Artifact Requirements](#artifact-requirements)
  - [Resource requirements](#resource-requirements)
  - [Evaluation runtime](#evaluation-runtime)
- [Structure and Content](#structure-and-content)
- [Getting Started](#getting-started)
  - [Load OVA](#load-ova)
  - [Replicating Experiments](#replicating-experiments)
    - [Examples Run by BdryReach](#examples-run-by-bdryreach)
    - [Examples Run by CORA](#examples-run-by-cora)
- [Functional badge](#functional-badge)
- [Reusable badge](#reusable-badge)
  - [License](#license)
  - [Tool Installation in a New Computer](#tool-installation-in-a-new-computer)
    - [Eigen3](#eigen3)
    - [Python](#python)
    - [boost](#boost)
    - [glpk](#glpk)
    - [Capd](#capd)
    - [BdryReach Toolkit Installation and Compilation of Test Cases](#bdryreach-toolkit-installation-and-compilation-of-test-cases)
  - [Simple user Guide](#simple-user-guide)
    - [Include Files](#include-files)
    - [Definition of Differential Equations](#definition-of-differential-equations)
    - [Parameter Configuration for Computing Reachable Sets](#parameter-configuration-for-computing-reachable-sets)
    - [Invoking the Boundary-based Method for Computing the Inner-approximations of Reachable Sets](#invoking-the-boundary-based-method-for-computing-the-inner-approximations-of-reachable-sets)
      - [Inner-approximation of Reachable Set Computation Interface](#inner-approximation-of-reachable-set-computation-interface)
    - [The Plotting of Results](#the-plotting-of-results)
    - [Compile and Run](#compile-and-run)
    - [Results Display](#results-display)
  
#  Getting Started

We provide the Open Virtual Appliance (OVA) image which can quickly use our tool without any additional installation and configuration. The virtual machine is installed with Ubuntu 20.04.6 LTS, **the username is laode and the password is 123456.**
#  Load OVA
Open the ova file **Ubuntu_aeOS.ova** in VMware Workstation (or other virtual machines). The codes and executable files are located in the director **./home**. The source code of BdryReach is located in **./BdryReach_code** and the source code of CORA is located in **./CORA_code**. (Herein we see **./home** as the root directory.)
##  Replicating Experiments 
**The examples files are located in ./examples, which contains all the examples showcased in Table 3-5. It has tow folders, ./examples/BdryReach and ./examples/CORA, which include examples run by BdryReach (our tool) and CORA (compared tool) respectively.**

### Examples Run by BdryReach
For example, if you want to replicate the examples in Table 3 run by our tool, then you can run **runTable3_Examples.sh** in the directory **./examples/BdryReach**. 
``` bash
./runTable3_Examples.sh
```

The output will be saved in the file **Table3_examples_output.txt** in the directory **./examples/BdryReach/Table3_examples** and the figures are stored in the directory **./examples/BdryReach/Table3_examples/figure**. Here is the some part output of **runTable3_Examples.sh**. It will output 
 the zonotopic inner-approximation (represented by a matrix, the first column denotes the center of the zonotope and the other columns denote its generators) computed by our tool and the computation time for each example.
```
Table3_examples

ElectroOsc

zonotopic inner-approximation:
[-6.678246781374445 0.07556403903050581 3.861000979026132e-05 -0.2942877976311855 ;
2.822964194048449 0.007304956271076335 0 0.07559828881651579 ;]
time cost: 22.32342

......
```

If you want to recompile the executable files in **./examples/BdryReach**, you should run the following commands in the directory **./BdryReach_code**.

``` bash
mkdir build
cd build
cmake ..
make
```
All the executable files will be built into the directories **./BdryReach_code/build/Table3_examples**, **./BdryReach_code/build/Table4_examples**, **./BdryReach_code/build/Table5_examples_shortTime**, **./BdryReach_code/build/Table5_examples_longTime**.

### Examples Run by CORA
If you want to replicate the examples in Table 3 run by the compared tool CORA and compute $\gamma_{min}$ of inner-approximations computed by BdryReah and CORA (since CORA cannot compute inner-approximations of examples in Table 4 and examples with long time horizon in Table 5, we only compute $\gamma_{min}$ of results computed by our tool in their corresponing scripts), then you can open matlab in the terminal,

```bash
/usr/local/Polyspace/R2020a/bin/matlab
```

next open and run the file **runTable3_Examples.m** in the directory **./examples/CORA**. this program will output 
* the computation time of each example run by CORA;
* the $\gamma_{min}$ of inner-approximation of each example computed by our tool BdryReach;
* the $\gamma_{min}$ of inner-approximation of each example computed by CORA.

The output will be saved in the file **Table3_examples_output.txt** in the directory **./examples/CORA/Table3_examples** and the figures are stored in the directory **./examples/CORA/Table3_examples/figure**. Here is the output of **runTable3_Examples.m**. 
```
Table3_examples

ExampleElectroOsc
Computation time: 41.7216 seconds  // the computation time run by CORA
Precision_BdR: 0.87911      // the precision of inner-approximation computed by our tool BdryReach
Precision_CORA: 0.57703     // the precision of inner-approximation computed by CORA

ExampleRoessler
Computation time: 42.5504 seconds
Precision_BdR: 0.76085
Precision_CORA: 0.78668

ExampleLotkaVolterra
Computation time: 288.2062 seconds
Precision_BdR: 0.6499
Precision_CORA: 0.34126

ExampleTank6
Computation time: 254.2386 seconds
Precision_BdR: 0.81845
Precision_CORA: 0.63015

ExampleBiologicalSystemI
Computation time: 126.8434 seconds
Precision_BdR: 0.9635
Precision_CORA: 0.89848

ExampleBiologicalSystemII
Computation time: 254.6293 seconds
Precision_BdR: 0.93196
Precision_CORA: 0.87238

ExampleTank12
Computation time: 2131.1678 seconds
Precision_BdR: 0.77182
Precision_CORA: 0.55813
```

# Functional badge

Since the $\gamma_{min}$ is computed by simulation method, the reuslts should be correct if the $\gamma_{min}$ is around the corresponing $\gamma_{min}$ recorded in the table ($\pm 0.05$). As for the computation time, it should be correct if it is proportional to the corresponding time recorded in the table.  

# Reusable badge

## License
This project is licensed under the GNU GPLv3 License - see the [LICENSE](./ReadMe%20file/LICENSE.md) file for details.

## Tool Installation in a New Computer
**We also provide the instruction on how to deploy our tool in a new computer with Linux OS.**

### Cmake
```bash
sudo apt-get install cmake
```

###  Eigen3

```bash
sudo apt-get install libeigen3-dev
```

### pip
```bash
sudo apt-get install pip
```

### Numpy

```bash
sudo apt-get install python-numpy
```

### Matplotlib

```
pip install matplotlib
```

###  boost

```bash
sudo apt-get install libboost-all-dev
```

###  glpk

```bash
sudo apt-get install libglpk-dev
```
### c\c++ compiler (optional)

```bash
sudo apt install build-essential
```
###  Capd
In the directory **./capd**.
```bash
sudo apt install libtool
autoreconf --install
cd ..
mkdir capd_nogui
cd capd_nogui
../capd/configure --prefix=/usr/local --without-gui --without-mpfr
make
sudo make install
```


### BdryReach Toolkit Installation and Compilation of Test Cases

in the directory **./BdryReach_code/**.

```bash
mkdir build
cd build
cmake ..
make
```



## Simple user Guide 

In this section we will show how to compute an inner-approximation using our tool BdryReach. The main source code is located in ./BdryReach_code/include/underApprox/underApprox.h. We take the example located in ./BdryReach_code/Table3_examples/
ElectroOsc.cpp as an instance.

### Include Files

```cpp
#include <plotter/matplotlibcpp.h>   // Header for computing reachable set outer-approximation
#include <plotter/plotter.h>          // Header for result visualization
#include <underApprox/underApprox.h>  // Header for includes the interface for computing reachable sets under approximation.
```

### Definition of Differential Equations

**We use the Capd library to define the form of the differential equations. Refer to the Capd documentation on [differential equation systems](https://capd.sourceforge.net/capdDynSys/docs/html/maps.html). Notably, the computation of our method requires validation of the obtained reachable set inner-approximation. Therefore, an additional definition for a time-inverted differential equation is necessary.**

```cpp
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -in[1];
    out[1] = -(0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

int dimIn = 3;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
```

### Parameter Configuration for Computing Reachable Sets

**We adopt parameter definitions similar to the MATLAB Reachability Analysis Toolbox CORA. For detailed meanings, refer to CORA's documentation.**
```cpp
NonlinearSys<double> mysys(f, 2, 0, 2);
NonlinearSys<double> mysysBack(fBack, 2, 0, 2);

ReachOptions<double> options;

// create R0
Vector_t<double> center(2);
center << 0, 3;
// center << -7.5, 2.8;
Matrix_t<double> generators(2,2);
generators<< 0.1,0,
                0,0.1;

Zonotope<double> R0_(center,generators);

options.set_time_step(0.2);
options.set_taylor_terms(4);
options.set_zonotope_order(50);
options.set_intermediate_order(50);
options.set_error_order(20);
options.set_alg("lin");
options.set_tensor_order(3);

options.set_tFinal(3);
options.set_tStart(0);

options.set_R0(R0_);

options.set_usekrylovError(1);
options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
```

### Invoking the Boundary-based Method for Computing the Inner-approximations of Reachable Sets

**This step invokes our boundary-based method for computing inner-approximations of reachable sets. Please refer to the interface for the meanings of various parameters.**
```cpp
vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 3, 1, 0.01, 0.005, 0.005, 50,50);
```
#### Inner-approximation of Reachable Set Computation Interface 
```cpp
template <typename Number>
        static vector<Zonotope<Number>> underReachClp(NonlinearSys<Number> mysys, 			
        NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double 	
        overRtime, int steps, double radius, double over_step, double bound_step, int Zover_order)
```
**Parameters:**
* **mysys:** ordinary differential equation for computing reachable sets.
* **mysysBack:** time-inverted ordinary differential equation for result verification.
* **options:** relevant configuration for outer-approximation of reachable set computation in the program.
* **R0:** initial set.
* **overRtime:** step size for inner-approximation of reachable set computation at each step.
* **steps:** number of iterations for inner-approximation of reachable set computation.
* **radius:** maximum allowed length of generator for facets.
* **over_step:** step size for outer-approximation computation for the entire set at each step in inner-approximation computation.
* **bound_step:** step size for outer-approximation computation for the boundary of the set at each step in inner-approximation computation.
* **Zover_order:** limit on the zonotope order for outer-approximation computation for the entire set at each step in inner-approximation computation.
  
### The Plotting of Results

For plotting the graphical results, we utilize the lightweight plotting library **Matplotlib for C++**." For specific usage instructions,please refer to [Matplotlib for C++ Documentation](https://matplotlib-cpp.readthedocs.io/en/latest/index.html).
```cpp
plt::figure_size(1200, 780);
for(int i = 1; i < underR.size(); i++){
    Plotter::plotZonotope(underR[i], 1, 2, "g");
}
Plotter::plotZonotope(R0_, 1, 2, "k");
plt::show();
```

### Compile and Run 

We use **cmake** to compile the program. 

Firstly, in the file **./BdryReach_code/CMakeLists.txt**, add the code ```add_subdirectory(Table3_examples) ```
to the last row.

* **./BdryReach_code/CMakeLists.txt** 
```cmake
cmake_minimum_required(VERSION 3.14)
project(reachSolver)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)

add_subdirectory(Table3_examples)

```
Next, in the file **./BdryReach_code/Table3_examples/CMakeLists.txt**, set the ```TARGET_NAME``` as the cpp file name ```ElectroOsc```, i.e., the code ```set(TARGET_NAME ElectroOsc)```.  
* **./BdryReach_code/Table3_examples/CMakeLists.txt** 
```cmake
cmake_minimum_required(VERSION 3.14)
project(reach_solver_test)

set(CMAKE_CXX_STANDARD 14)

set(TARGET_NAME ElectroOsc) # Set the target name
#ElectroOsc
#Roessler
#LotkaVolterra
#Tank6
#BiologicalSystemI
#BiologicalSystemII
#Tank12
find_package(Eigen3 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED QUIET)

add_executable(${TARGET_NAME} ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} PUBLIC
        ${PYTHON_LIBRARIES}
        /usr/local/lib/libcapd.so
        Eigen3::Eigen
        glpk
        reach_solver)
```

in the directory **./BdryReach_code**, run the following bash code, 
```bash
mkdir build
cd build
cmake ..
make
```
if the compile successes, you can find an executable file in the directory **/BdryReach_code/build/Table3_examples**, then you can run the executable file **ElectroOsc**.

```bash
cd Table3_examples
./ElectroOsc
```
### Results Display

The output will show the zonotopic inner-approximation (represented by a matrix, the first column denotes the center of the zonotope and the other columns denote its generators) computed by our tool. Additional, the output will show the computation time of the whole program. 
```
zonotopic inner-approximation: //the computed inner-approximation 
[-6.678246781374445 0.07556403903050581 3.861000979026132e-05 -0.2942877976311855 ;
2.822964194048449 0.007304956271076335 0 0.07559828881651579 ;]
time cost: 22.754546   
```
The output will also show the images of computed zonotopes at each iteration on $x_1$- $x_2$ plane. 
<p align="center">
  <img src=./ReadMe%20file/BdryReach_ElectorOsc.svg>
</p>

The black region represents the initial set , the green region represents the inner-approximation of the reachable set, while, for comparison, the blue region represents the outer-approximation of the reachable set.

