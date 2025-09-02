# FluidX3D-Cd-tool
Drag coefficients (C_d) visualization and analysis tool for FluidX3D by Dr. Moritz Lehmann
https://github.com/ProjectPhysX/FluidX3D 


<p align="center">
  <img src="misc/Re10000_Q_10x.gif" alt="Simulation demo (10 times faster" width="45%"/>
  <img src="misc/Re10000_Cd_10x.gif" alt="Cd output demo (10 times faster)" width="45%"/>
</p>

<!-- 
![Simulation demo (10 times faster)](misc/Re_10000_Q_10x.gif)
![Cd output demo (10 times faster)](misc/Re_1000_Cd_10x.gif)
--!>

## Guideline

copy live_cd_plot.py to the directory where FluidX3D.exe is located (\FluidX3D-master\bin)
that is also where the .dat file is generated in its subdirectory (eg. FP16C\500MB\Cd.dat)

## STL files

sphere_1000mm.stl
MosquitoSolo.stl
rambling_wreck_scaled_ASCII.stl

*note that rambling wrek needs a different BC setting by adding wall (ground) in the bottom. 

**set:** 

