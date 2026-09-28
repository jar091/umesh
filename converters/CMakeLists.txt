# ======================================================================== #
# Copyright 2018-2022 Ingo Wald                                            #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
#     http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
# ======================================================================== #

add_subdirectory(vtu)
add_subdirectory(bhac2umesh)
add_subdirectory(foam2umesh)

# ------------------------------------------------------------------
# tools: big lander (the large nasa lander that comes in multiple
# individual files)
# ------------------------------------------------------------------
add_executable(umeshImportLanderFun3D
  importLanderFun3D.cpp
  )
target_link_libraries(umeshImportLanderFun3D
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# tools: read older umesh variant with other magic number, and write
# as a current-version umesh
# ------------------------------------------------------------------
# add_executable(umeshImportOldUMesh
#   importOldUMesh.cpp
#   )
# target_link_libraries(umeshImportOldUMesh
#   PUBLIC
#   umesh
#   )

# ------------------------------------------------------------------
# tool that reads a 'raw' structured volume file, chops it into little
# grids, and creates a 'umesh' that has these grids
# ------------------------------------------------------------------
add_executable(umeshRawToGrids
  rawToGrids.cpp
  )
target_link_libraries(umeshRawToGrids
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# tool that takes a fun3d model (with mesh and volume files) and
# converts them to one umesh per rank and one raw scalars per rank per
# variable per time step
# ------------------------------------------------------------------
add_executable(umeshBreakApartFun3D
  breakApartFun3D.cpp
  )
target_link_libraries(umeshBreakApartFun3D
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# tool that takes a fun3D input that comes in the form of
#
# a) a single ".lb8" ugrid64 file for the mesh, and
#
# b) a set of per-rank scalar files ...
#
# and creates a single umesh with a single variable
# ------------------------------------------------------------------
add_executable(fun3DToUmesh
  fun3DToUmesh.cpp
  )
target_link_libraries(fun3DToUmesh
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# tool that takes a fun3D input that comes in the form of
#
# a) a single ".umesh" file for the mesh, typically one of the
#    per rank .umesh files that resulted from running a partitioner
#
# b) a set of per-rank scalar files ...
#
# and creates a single umesh with a single variable
# ------------------------------------------------------------------
add_executable(fun3DExtractVariable
  fun3DExtractVariable.cpp
  )
target_link_libraries(fun3DExtractVariable
  PUBLIC
  umesh
  )


