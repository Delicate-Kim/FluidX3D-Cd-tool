# LBM_lid_driven_flow_in_polygons
Lattice Boltzmann MATLAB code for 2D lid-driven flow in polygonal and complex cavities    
Author: *Soohwan Kim*  

---

## Method Overview
This project implements a **D2Q9 Lattice Boltzmann solver** with the **SRT (Single Relaxation Time)** model and the **BGK (Bhatnagar–Gross–Krook)** collision operator. The implementation closely follows *Mohamad (2011)* and incorporates a clearer distinction between physical, nondimensional, and numerical/discrete parameters, as discussed in *Latt (2008)*.

**D2Q9 velocity set:**
<pre> 
 c7  c3  c6 
