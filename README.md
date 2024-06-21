<div align="center">
  	<h1>
    	BdryReach User Guide
  	</h1>
	<br />
    <div>
		<a href="https://github.com/BdryReach/BdryReach">
			<img width="428" src="result_picture/BdryReach.svg" alt="BdryReach System">
		</a>
	</div>
</div>

## Introduction
 Reachability analysis, a pivotal component in formal verification of dynamic systems, entails the computation of outer- and inner-approximations of reachable sets. Outer-approximations, supersets of the exact reachable sets, are instrumental in safety verification, as a system is deemed safe if the computed outer-approximation excludes unsafe states. Conversely, inner-approximations, subsets of the exact reachable sets, confirm the existence of trajectories that intersect with unsafe states. The propagation of the initial set often introduces computation error accumulation, a.k.a. wrapping effect, which greatly affects the computation accuracy of the outer-/inner- approximations. To address this, Xue et al. introduced a novel <a href="http://lcs.ios.ac.cn/~xuebai/publication.html"><strong>set-boundary propagation based reachability method</strong></a>, leveraging a deep investigation on the topological structure. This approach surpasses the common practice of partitioning the initial set, as it reduces computational burden by focusing on the boundary of the initial set rather than its entirety. This tool implements the set-boundary reachability method, in which the computed outer- and inner-approximations are represented by zonotopes.

## 1. Installation
To run BdryReach in a Linux system, it is required to install the **cmake** tool and the following third-party libraries.

### 1.1 Required Third-Party Libraries

| Library | Website | Version |
| --- | --- | --- |
| Eigen | [http://eigen.tuxfamily.org/index.php?title=Main Page](http://eigen.tuxfamily.org/index.php?title=Main%20Page) | 3.34 |
| Python | [https://www.python.org/](https://www.python.org/) | 3.8 |
| Capd | [http://capd.ii.uj.edu.pl/](http://capd.ii.uj.edu.pl/) | latest |
| Boost | [https://www.boost.org/](https://www.boost.org/) | 1.67.0.0 |
| GLPK | [https://www.gnu.org/software/glpk/](https://www.gnu.org/software/glpk/) | 4.65-2 |
| Pip | [https://pypi.org/project/pip/](https://pypi.org/project/pip/) | latest |
| Numpy | [https://numpy.org/](https://numpy.org/) | latest |
| Matplotlib | [https://matplotlib.org/](https://matplotlib.org/) | latest |
### 1.2 Installation for Dependencies

#### 1.2.1 Eigen
```bash
sudo apt-get install libeigen3-dev
```
#### 1.2.2 Pip
```bash
sudo apt-get install pip
```
#### 1.2.3 Numpy

```bash
sudo apt-get install python-numpy
```

#### 1.2.4 Matplotlib

```
pip install matplotlib
```

#### 1.2.5  Boost

```bash
sudo apt-get install libboost-all-dev
```

#### 1.2.6  GLPK

```bash
sudo apt-get install libglpk-dev
```
#### 1.2.7 Capd
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
### 1.3 Download of BdryReach and Compilation of Test Cases
```bash
git clone https://github.com/ASAG-ISCAS/BdryReach.git
cd BdryReach/
mkdir build
cd build
cmake ..
make
```
## 2. Usage

### 2.1 Interface for Computing Outer- and Inner-approximations 
### 2.1.1  Interface for Computing Outer-approximations
```cpp
template <typename Number>
static vector<ReachableSet<Number>> BdReach(NonlinearSys<Number> mysys, ReachOptions<Number> options, Zonotope<Number> R0)
```
**Parameters:**
* **mysys:** differential equations.
* **options:** configuration for outer-approximation of reachable set computation.
* **R0:** initial set.


### 2.1.2 Interface for Computing Inner-approximations
```cpp
template <typename Number>
        static vector<Zonotope<Number>> underReachClp(NonlinearSys<Number> mysys, 			
        NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double 	
        overRtime, int steps, double radius, double over_step, double bound_step, int Zover_order)
```
**Parameters:**
* **mysys:** differential equatiions.
* **mysysBack:** time-inverted differential equations.
* **options:** configuration for the computation of outer-approximations in the program.
* **R0:** initial set.
* **overRtime:** time step size for the comoutation of inner-approximations.
* **steps:** number of iterations for for the comoutation of inner-approximations.
* **radius:** maximum allowed length of generators of facets.
* **over_step:** time step size for computing the outer-approximation of the entire reachable set at each step.
* **bound_step:** time step size for computing the outer-approximation of the boundary of reachable set at each step
* **Zover_order:** limit on the zonotope order in the computation of the outer-approximation of the entire reachable set at each step.
### 2.2 Example for Computing Outer-approximations
**As an example, we perform the computation of outer-approximations for the VanderPol model. The file computes the outer-approximations from the initial region ([1.23, 1.57], [2.34, 2.46]) over the time interval [0, 6.74] (in seconds).The specific file location is:**
```RobotFramework
/examples/overVanderPol.cpp.
```
### 2.2.1 Header Files
```cpp
#include <overApprox/overApprox.h> // Header file with interfaces for computing outer-approximations.
#include <plotter/matplotlibcpp.h> // Header file for Matplotlib C++ plotting library.
#include <plotter/plotter.h> // Header file for result plotting.
```
### 2.2.2 Definition of Differential Equations
**We define the form of differential equations using the Capd library. For detailed information on the differential equations in Capd, please refer to the [Capd documentation](https://capd.sourceforge.net/capdDynSys/docs/html/maps.html).**


```cpp
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0]/*+ in[2]*/;
}

// Input dimension of the differential equations
int dimIn = 3; // The input dimension of the differential equations. Since the default input includes control u, the input dimension is one more than the output dimension.

// Output dimension of the differential equations
int dimOut = 2; // The output dimension of the differential equations.

// Parameter settings for the differential equations. Since this differential equations have no parameters, it is set to 0.
int noParam = 0;

// Maximum order for Taylor expansion of the differential equations.
int MaxDerivativeOrder = 3; // The maximum order to which the differential equations is expanded using Taylor series.

// Creating IMap for interval computations
IMap f(_f, dimIn, dimOut, noParam, MaxDerivativeOrder); // Constructing IMap for interval Computations.


```
### 2.2.3 Parameter Configuration for Computing Outer-approximations
**Here, we adopt the same parameter settings as the Continuous Reachability Analyzer CORA. The meanings of each parameter can be found in CORA's documentation. Please refer to the [manual of CORA](result_picture/Cora2021Manual.pdf).**
```cpp
    NonlinearSys<double> mysys(f, 2, 0, 2);
    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(2);
    Vector_t<double> c(2);
    c << 1.4, 2.4;
    Matrix_t<double> generators(2,1);
    Matrix_t<double> G(2,2);
    G<< 0.17,0,
                 0,0.06;
    Zonotope<double> R0_(c,G);

    center << 1.4, 2.46;
    generators<< 0.17,
                 0;

    options.set_R0(R0_);

    options.set_time_step(0.005);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tFinal(6.74);
    options.set_tStart(0);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
```
### 2.2.4 Invocation of the Boundary-based Method for Computing Outer-approximations
This step invokes our boundary-based method for computing outer-approximations of reachable sets. Please refer to **Section 2.1.1** for the meanings of  parameters.
```cpp
vector<ReachableSet<double>> BdReachset = OverApprox::BdReach(mysys, options, R0_);
```
### 2.2.5 The Plotting Results
For plotting the results, we utilize the lightweight plotting library **Matplotlib for C++**." For specific usage instructions, please refer to [Matplotlib for C++ Documentation](https://matplotlib-cpp.readthedocs.io/en/latest/index.html).
```cpp
plt::figure_size(1200, 780);
for(int i = 0; i < BdReachset.size(); i++){
    Plotter::plotReach(BdReachset[i], 1, 2, "b");
}
plt::show();
```
### 2.2.6 Results Display
**We employ both BdryReach and CORA to compute the outer-approximations of the reachable set starting from the initial region ([1.23, 1.57], [2.34, 2.46]) over the time interval [0, 6.74] (in seconds). The blue region represents the results obtained by BdryReach, while the red region corresponds to the results of CORA. It is evident that the outer-approximations computed by BdryReach exhibit significantly higher accuracy compared to CORA.**
<p align="center">
  <img src=result_picture/2.2.6.png>
</p>

## 2.3 Test Case for Inner-approximating Reachable sets
**We also take the inner-approximation of the reachable set computation for the VanderPol model as an example. The file computes the inner-approximation from the initial region ([1.23, 1.57], [2.34, 2.46]) with a step size of 0.1s over the time interval [0, 0.8] (in seconds). The specific file location is /examples/underVanderPol.cpp.**
### 2.3.1 Include Files

```cpp
#include <plotter/matplotlibcpp.h>   // Header for computing reachable set outer-approximation
#include <plotter/plotter.h>          // Header for result visualization
#include <underApprox/underApprox.h>  // Header for includes the interface for computing reachable sets under approximation.
```
### 2.3.2 Defining Differential Equations

**We use the Capd library to define the form of the differential equations (refer to the Capd documentation on [differential equation systems](https://capd.sourceforge.net/capdDynSys/docs/html/maps.html)). Notably, the computation of our tool requires validation of the obtained reachable set inner-approximation. Therefore, an additional definition for a time-inverted differential equation is necessary.**

```cpp
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    // Forward differential equation
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0];
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    // Reverse differential equation
    out[0] = -in[1];
    out[1] = -(mu * (1-in[0]*in[0])*in[1] - in[0]);
}

int dimIn = 3;
int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn, dimOut, noParam, MaxDerivativeOrder);
IMap fBack(_fBack, dimIn, dimOut, noParam, MaxDerivativeOrder);
```
### 2.3.3 Parameter Configuration for Computing Reachable Sets

**We adopt parameter definitions similar to CORA. For detailed meanings, refer to CORA's documentation.**
```cpp
NonlinearSys<double> mysys(f, 2, 0, 2);
NonlinearSys<double> mysysBack(fBack, 2, 0, 2);

ReachOptions<double> options;

// Create R0
Vector_t<double> center(2);
center << 1.4, 2.4;
Matrix_t<double> generators(2,2);
generators << 0.17, 0,
              0, 0.06;

Zonotope<double> R0_(center, generators);

options.set_taylor_terms(4);
options.set_zonotope_order(50);
options.set_intermediate_order(50);
options.set_error_order(20);
options.set_alg("lin");
options.set_tensor_order(3);

options.set_tFinal(6);
options.set_tStart(0);

options.set_R0(R0_);

options.set_usekrylovError(1);
options.set_max_error(DBL_MAX * Eigen::MatrixXd::Ones(2,1));
```
### 2.3.4 Invocation of the Boundary-based Method for Computing the Inner-approximations of Reachable Sets

**This step invokes our boundary-based method for computing inner-approximations of reachable sets. Please refer to Section 2.1.1 for the meanings of various parameters.**
```cpp
vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 0.1, 8, 0.01, 0.05, 0.01, 50);
```
### 2.3.5 Results Plotting

For plotting the graphical results, we utilize the lightweight plotting library **Matplotlib for C++**." For specific usage instructions,please refer to [Matplotlib for C++ Documentation](https://matplotlib-cpp.readthedocs.io/en/latest/index.html).
```cpp
plt::figure_size(1200, 780);
for(int i = 1; i < underR.size(); i++){
    Plotter::plotZonotope(underR[i], 1, 2, "g");
}
Plotter::plotZonotope(R0_, 1, 2, "k");
plt::show();
```
### 2.3.6 Results Display

**We computed the inner-approximation of the reachable set for the VanderPol model starting from the initial region ([1.23, 1.57], [2.34, 2.46]) with a step size of 0.1s over the time interval [0, 0.8] (in seconds). The green region represents the inner-approximation of the reachable set, while, for comparison, the blue region represents the outer-approximation of the reachable set. It is evident that the inner-approximation computed provides tight results.**
<p align="center">
  <img src=result_picture/2.3.6.png>
</p>

## 3 Reproducing Results

**All experiments were conducted on a virtual machine system with an Intel® Core™ i7-10750H CPU @ 2.60GHz × 8, and 5.8 GiB of memory, running Ubuntu 20.04.3 LTS. The experimental files are located in the /examples directory. Please note that different runtime environments may introduce some variations in the results. We provide 8 C++ files to replicate all experiments. To run these experiments, simply execute the following CMake statement in the /BdaryReach directory.**
```bash
mkdir build
cd build
cmake ..
make

./examples/underBiological
./examples/underTank
./examples/underVanderPol
./examples/underElectroOsc
./examples/overElectroOsc
./examples/overTank
./examples/overVanderPol
./examples/over2Rrobot
```
# Acknowledgement

This project was completed by Dejin Ren, Changyuan Zhao, and Bai Xue from the Institute of Software, Chinese Academy of Sciences (ISCAS).
