## 3D Graphics Rasterizer  
Rasterization is a rendering algorithm for computer graphics. It takes a 3D representation of geometry and maps it to 2D space to be rendered on the screen. 

## Running it  
To compile, you must have cmake and vtk installed.  
- cmake, at least version 3.6  
- vtk, version 6.\*  

To execute, clone the repo and run the following commands  
`cmake .`  
`make`  
Navigate the directory created `rasterizer.app` to find the executable. Execute with  
`./rasterizer.app/.../rasterizer filename`  
where 'filename' is the geometry file.  
The geometry file must be a '.vtk' file in ASCII format with 'PolyData' representation.  
## Caveats  
- Compilation has been tested on MacOS 10.11; cmake offers cross-platform compilation, but this is untested.  
