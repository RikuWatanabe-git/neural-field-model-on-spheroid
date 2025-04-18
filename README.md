## **Neural Field Model on Oblate Spheroids**

Solve neural field equations on oblate spheroids.

### **Description**

This repository provides a numerical solver for neural field equations on oblate spheroids with small flattening.  
For details on the mathematical formulation and implementation, see Ref. [1].

The code is based on the open code in Ref.[2].

The surface discretization and area element computations on the triangular mesh are based on the method proposed in Ref. [3].

### **Requirements**

#### **Julia**
- PyCall
- PyPlot
- DelimitedFiles
- LinearAlgebra
- Colors

#### **Python**
- NumPy  
- SciPy
- geographiclib  

Ensure these dependencies are installed before running the scripts.

#### **Usage**

First, run distance_calculator.jl, then run spheroid.jl. Model parameters can be changed by editing these scripts.

## References
1. H. Ishii and R. Watanabe <br>
"Spot solutions to a neural field equation on oblate spheroids" (preprint) <br>

2. R. Nishide and S. Ishihara <br>
"Pattern Propagation Driven by Surface Curvature"<br>
Phys. Rev. Lett. 128, (2022), 224101 <br>

3. G. Xu <br>
"Discrete Laplace-Beltrami operator on sphere and optimal spherical triangulations"<br>
Int. J. Comput. Geom. Appl. 16, pp.75-93 (2006) <br>
