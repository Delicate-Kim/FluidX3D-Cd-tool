# FluidX3D-Cd-tool
Drag coefficients (C_d) visualization and analysis tool for FluidX3D by Dr. Moritz Lehmann
https://github.com/ProjectPhysX/FluidX3D 

This repository provides helper scripts and example to output drag coefficients and visually track the values over timesteps. It also provides the baseline measurement for drag coefficients of a sphere.

Example of sphere 

<p align="center">
  <img src="misc/Re_10000_Q_10x.gif" alt="Simulation demo (10 times faster" width="45%"/>
  <img src="misc/Re10000_Cd_10x.gif" alt="Cd output demo (10 times faster)" width="45%"/>
</p>

<!-- 
![Simulation demo (10 times faster)](misc/Re_10000_Q_10x.gif)
![Cd output demo (10 times faster)](misc/Re_10000_Cd_10x.gif)
--!>

## Guideline

copy live_cd_plot.py to the directory where FluidX3D.exe is located (\FluidX3D-master\bin). This is also where the .dat file is generated in its subdirectory (eg. FP16C\500MB\Cd.dat)
copy the stl files where you call the stl files (typically, \FluidX3D-master\stl)
replace your "void main setup() {...}" with the snippet below to conduct sphere drag scneario.

[code format]
void main_setup() { // Ahmed body; required extensions in defines.hpp: FP16C, FORCE_FIELD, EQUILIBRIUM_BOUNDARIES, SUBGRID, optionally INTERACTIVE_GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	
	const uint memory = 500u; // available VRAM of GPU(s) in MB, you can increase even higher
	const float lbm_u = 0.05f;
	const float box_scale = 4.0f; //6.0f
	const float si_u = 0.148; // 0.000148 = Re 10
	const float si_nu = 1.48E-5f, si_rho = 1.225f;
	const float si_D = 1.0f;
	const float PI = 3.14159265358979323846;
	const float si_A = 0.25f * si_D * si_D * PI;
	const float si_T = 20000000000000.0f;
	const float si_Lx = box_scale * si_D;
	const float si_Ly = 2.0f * box_scale * si_D;
	const float si_Lz = box_scale * si_D;
	
	const uint3 lbm_N = resolution(float3(si_Lx, si_Ly, si_Lz), memory); // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
	units.set_m_kg_s((float)lbm_N.y, lbm_u, 1.0f, box_scale * si_D, si_u, si_rho);
	const float lbm_nu = units.nu(si_nu);
	const ulong lbm_T = units.t(si_T);
	const float lbm_length = units.x(si_D);
	print_info("Re = " + to_string(to_uint(units.si_Re(si_D, si_u, si_nu))));
	LBM lbm(lbm_N, lbm_nu);
	info.lbm = &lbm;
	// ###################################################################################### define geometry ######################################################################################
	Mesh* mesh = read_stl(get_exe_path() + "../stl/sphere_1000mm.stl", lbm.size(), lbm.center(), float3x3(float3(0, 0, 1), radians(90.0f)), lbm_length); //ahmed_25deg_m.stl
	mesh->translate(float3(0.0f, units.x(0.5f * (0.5f * box_scale * si_D - si_D)) - mesh->pmin.y, 0.0f));
	
	lbm.voxelize_mesh_on_device(mesh, TYPE_S | TYPE_X); // https://github.com/nathanrooy/ahmed-bluff-body-cfd/blob/master/geometry/ahmed_25deg_m.stl converted to binary
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);

	if (lbm.flags[n] != TYPE_S) lbm.u.y[n] = lbm_u;
	if (x == 0u || x == Nx - 1u || y == 0u || y == Ny - 1u || z == 0u || z == Nz - 1u) lbm.flags[n] = TYPE_E;
		}); // ####################################################################### run simulation, export images and data ##########################################################################
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE | VIS_FIELD;
	lbm.graphics.field_mode = 1;
	lbm.graphics.slice_mode = 1;
	//lbm.graphics.set_camera_centered(20.0f, 30.0f, 10.0f, 1.648722f);
	lbm.run(0u, lbm_T); // initialize simulation
#if defined(FP16S)
	const string path = get_exe_path() + "FP16S/" + to_string(memory) + "MB/";
#elif defined(FP16C)
	const string path = get_exe_path() + "FP16C/" + to_string(memory) + "MB/";
#else // FP32
	const string path = get_exe_path() + "FP32/" + to_string(memory) + "MB/";
#endif // FP32
	lbm.write_status(path);
	write_file(path + "Cd.dat", "# t\tCd\n");
	const float3 lbm_com = lbm.object_center_of_mass(TYPE_S | TYPE_X);
	print_info("com = " + to_string(lbm_com.x, 2u) + ", " + to_string(lbm_com.y, 2u) + ", " + to_string(lbm_com.z, 2u));
	while (lbm.get_t() <= lbm_T) { // main simulation loop
		Clock clock;
		const float3 lbm_force = lbm.object_force(TYPE_S | TYPE_X);
		const float Cd = units.si_F(lbm_force.y) / (0.5f * si_rho * sq(si_u) * si_A); // expect Cd to be too large by a factor 1.3-2.0x; need wall model
		print_info("Cd = " + to_string(Cd, 3u) + ", t = " + to_string(clock.stop(), 3u));
		write_line(path+"Cd.dat", to_string(lbm.get_t())+"\t"+to_string(Cd, 3u)+"\n");
		lbm.run(1u, lbm_T);
	}
	//lbm.write_status(path);
} /**/

## drag coefficient sweeping over Reynolds number

This is a result by testing different Reynolds number using memory of 500MB
Testing device:
| Device Name    | NVIDIA GeForce RTX 4070 Laptop GPU                         |
| Device Vendor  | NVIDIA Corporation                                         |
| Device Driver  | 566.07 (Windows)                                           |
| OpenCL Version | OpenCL C 3.0                                               |
| Compute Units  | 36 at 1980 MHz (4608 cores, 18.248 TFLOPs/s)               |
| Memory, Cache  | 8187 MB VRAM, 1008 KB global / 48 KB local                 |
| Buffer Limits  | 2046 MB global, 64 KB constant


Example variable for Re = 10,000
| Grid Resolution |                                 158 x 315 x 158 = 7863660 |
| Grid Domains    |                                             1 x 1 x 1 = 1 |
| LBM Type        |                                    D3Q19 SRT (FP32/FP16C) |
| Memory Usage    |                                 CPU 217 MB, GPU 1x 508 MB |
| Max Alloc Size  |                                                    284 MB |
| Time Steps      |                                          4662000168730624 |
| Kin. Viscosity  |                                                0.00039375 |
| Relaxation Time |                                                0.50118124 |
| Reynolds Number |                                               Re < 566285 |

![Re vs Cd comparison with experiment and simulation](misc/Re_Cd_overlap.png)

The measured Cd values in solid black dots are plotted while overalaped with plot from (https://kdusling.github.io/teaching/Applied-Fluids/ImageDisplay.html?src=DragSphere)
Drag coefficient of a smooth sphere as a function of Reynolds number. Open blue triangles are data from Maxworthy, 1965 (add link: https://doi.org/10.1017/S002211206500143X). Open blue squares are data from Roos and Willmarth, 1971 (https://doi.org/10.2514/3.6164). The remaining blue data points are from earlier experiments, Schlichting, 1979 (https://doi.org/10.1007/978-3-662-52919-5).

note
- the simulation data was not perfectly matched to experimental data but still is within order of magnitude of 1 (O(10)) while following general trendline
- no drag crisis was observed near Re = 10^5 - 10^6
- resolution study may be worth to try

## STL files

sphere_1000mm.stl: sphere with diameter of 1000mm
MosquitoSolo.stl: free STL model from Cults3D (cults3d.com/en/3d-model/art/mosquito-model)
Ford_GT_2017.stl: free STL model from Thingiverse (www.thingiverse.com/thing:2045466)
Mustang_2014.stl: free STL model from Cults3D (www.thingiverse.com/thing:4978646)

*note that the mosquito wing does not flap and the the outcome drag coefficient is assumming fixed body.
*note that cars need a different BC setting by adding wall (ground) in the bottom. 


