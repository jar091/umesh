# ======================================================================== #
# Copyright 2018-2020 Ingo Wald                                            #
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

add_definitions(${OWL_CXX_FLAGS})

# ------------------------------------------------------------------
# print basic info of a umesh: #verts, tris, quads, tets, etc, as well
# as bounding box.
# ------------------------------------------------------------------
add_executable(umeshInfo
  info.cpp
  )
target_link_libraries(umeshInfo
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# merges several umesh files (e.g., neighbouring parts of a
# data-parallel decomposition) into one; shared vertices are not merged
# ------------------------------------------------------------------
add_executable(umeshMergeParts
  mergeParts.cpp
  )
target_link_libraries(umeshMergeParts
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# print basic info of a umesh: #verts, tris, quads, tets, etc, as well
# as bounding box.
# ------------------------------------------------------------------
add_executable(umeshReplaceScalarField
  replaceScalarField.cpp
  )
target_link_libraries(umeshReplaceScalarField
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# same as info, with added sanity checks (valid indices etc)
# ------------------------------------------------------------------
add_executable(umeshSanityCheck
  sanityCheck.cpp
  )
target_link_libraries(umeshSanityCheck
  umesh
  )
# ------------------------------------------------------------------
# import a OFF file to umesh format
# ------------------------------------------------------------------
add_executable(umeshImportOFF
  importOFF.cpp
  )
target_link_libraries(umeshImportOFF
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# load umesh(es) and save list of all bounding boxes
# ------------------------------------------------------------------
add_executable(umeshExportBBs
  exportBBs.cpp
  )
target_link_libraries(umeshExportBBs
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# import a RAW *structured* mesh as if it was a umesh with hexes
# ------------------------------------------------------------------
add_executable(umeshImportStructured
  importStructured.cpp
  )
target_link_libraries(umeshImportStructured
  PUBLIC
  umesh
  )

  # ------------------------------------------------------------------
# import a OFF file to umesh format
# ------------------------------------------------------------------
add_executable(umeshImportUGrid64
  importUGrid64.cpp
  )
target_link_libraries(umeshImportUGrid64
  PUBLIC
  umesh
  )
add_executable(umeshImportUGrid32
  importUGrid32.cpp
  )
target_link_libraries(umeshImportUGrid32
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# tetrahedralize a (general) umesh into a tet-only umesh, by
# introducing new vertices into quad centers
# ------------------------------------------------------------------
add_executable(umeshTetrahedralize
  tetrahedralize.cpp
  )
target_link_libraries(umeshTetrahedralize
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# given a umesh with volumetric elemnets (any sort), compute a new
# umesh (containing only triangles) that contains the triangular
# iso-surface for given iso-value. Input *must* have a per-vertex
# scalar field, but can have any combinatoin of volumetric
# elemnets; tris and quads in the input get ignored; input remains
# unchanged. 
# ------------------------------------------------------------------
add_executable(umeshExtractIsoSurface
  extractIsoSurface.cpp
  )
target_link_libraries(umeshExtractIsoSurface
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# computes the connectivity (tet and facets per face, and faces per tet) for a given tet mesh. only allowed for tet meshes
# ------------------------------------------------------------------
add_executable(umeshComputeTetConnectivity
  computeTetConnectivity.cpp
  )
target_link_libraries(umeshComputeTetConnectivity
  PUBLIC
  umesh
  )


# ------------------------------------------------------------------
# computes the outer shell of a tet-mesh, ie, all the triangle and/or
# bilinear faces that are _not_ shared between two neighboring
# elements
# ------------------------------------------------------------------
add_executable(umeshComputeShell
  extractShell.cpp
# do not use this any more:  computeShell.cpp
  )
target_link_libraries(umeshComputeShell
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# partitions a umesh, into either a BVH, or kd-tree. Note we don't
# actually save teh kdtree/bvh, only the final bricks, and the
# primrefs that refer into the original mesh.
# ------------------------------------------------------------------
add_executable(umeshPartitionSpatially
  partitionSpatially.cpp
  )
target_link_libraries(umeshPartitionSpatially
  PUBLIC
  umesh
  )

add_executable(umeshPartitionObjectSpace
  partitionObjectSpace.cpp
  )
target_link_libraries(umeshPartitionObjectSpace
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# no computations at all - just dumps the mesh that's implicit in the
# input umesh, and dumps it in obj format (other prims get ignored)
# ------------------------------------------------------------------
add_executable(umeshExtractSurfaceMesh
  extractSurfaceMesh.cpp
  )
target_link_libraries(umeshExtractSurfaceMesh
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# computes the outer shell faces of the volume elements (ie, all those
# faces that are NOT shared between two elements, then dumps those in
# obj or BTM format)
# ------------------------------------------------------------------
add_executable(umeshExtractShell
  extractShell.cpp
  )
target_link_libraries(umeshExtractShell
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# no computations at all - just dumps the mesh that's implicit in the
# input umesh, and dumps it in obj format (other prims get ignored)
# ------------------------------------------------------------------
add_executable(umeshExaToUMesh
  exaToUMesh.cpp
  )
target_link_libraries(umeshExaToUMesh
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# no computations at all - just dumps the mesh that's implicit in the
# input umesh, and dumps it in obj format (other prims get ignored)
# ------------------------------------------------------------------
add_executable(umeshFixNegativeVolumeElements
  fixNegativeVolumeElements.cpp
  )
target_link_libraries(umeshFixNegativeVolumeElements
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# computes faces, and face-connectiviy information, for a umesh; once
# via gpu, once via tbb
# ------------------------------------------------------------------
add_subdirectory(computeFaces)


# ------------------------------------------------------------------
# generates soe test cases for exa to umesh generator
# ------------------------------------------------------------------
add_executable(exaTestCases
  exaTestCase.cpp
  )
target_link_libraries(exaTestCases
  PUBLIC
  umesh
  )

# ------------------------------------------------------------------
# convert umesh to vtu
# ------------------------------------------------------------------
if (UMESH_USE_VTK)
  add_subdirectory(toVTU)
endif()
