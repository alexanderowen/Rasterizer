cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(rasterizer)
SET(VTK_DIR /Users/hank/Hartree/VTK/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(rasterizer MACOSX_BUNDLE rasterizer)

if(VTK_LIBRARIES)
  target_link_libraries(rasterizer ${VTK_LIBRARIES})
else()
  target_link_libraries(rasterizer vtkHybrid)
endif()

