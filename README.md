## **Neural Field Model on Oblate Spheroids**

Solve neural field equations on oblate spheroids.

### **Description**

This repository provides a numerical solver for the Amari neural field model on oblate spheroids with small flattening.  
For details on the mathematical formulation and implementation, see Ref. [1].

The code is based on the open code in Ref.[2].

The surface discretization and area element computations on the triangular mesh are based on the method proposed in Ref. [3].

### **Requirements**

#### **Julia (tested on Julia x.x.x)**
- PyPlot  
- DelimitedFiles  

#### **Python (tested on Python x.x.x)**
- NumPy  
- SciPy  
- Geographiclib  

Ensure these dependencies are installed before running the scripts.

#### **Usage**

Run distance_calculator.py first on Python, then run spheroid.jl on Julia.Model parameters can be changed by editing these scripts.

## Licence

[MIT](https://github.com/tcnksm/tool/blob/master/LICENCE)

## References
1. H. Ishii and R. Watanabe <br>
"Stationary spot solutions to a neural field equation on oblate spheroids" (preprint) <br>

2. R. Nishide and S. Ishihara <br>
"Pattern Propagation Driven by Surface Curvature"<br>
Phys. Rev. Lett. 128, (2022), 224101 <br>

3. G. Xu <br>
"Discrete Laplace-Beltrami operator on sphere and optimal spherical triangulations"<br>
Int. J. Comput. Geom. Appl. 16, p.75--93 (2006) <br>
