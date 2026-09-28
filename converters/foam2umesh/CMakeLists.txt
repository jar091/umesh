# ======================================================================== #
# Copyright 2018-2026 Ingo Wald                                            #
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

# ------------------------------------------------------------------
# decomposed OpenFOAM case (binary polyMesh + volScalar/VectorField)
# -> N stand-alone .umesh parts (see README.md)
# ------------------------------------------------------------------
find_package(OpenMP)
if (NOT OpenMP_CXX_FOUND)
  message("Skipping foam2umesh: OpenMP not found")
  return()
endif()

add_executable(foam2umesh
  main.cpp
  Convert.cpp
  Convert.h
  FoamIO.cpp
  FoamIO.h
  )
target_compile_features(foam2umesh PRIVATE cxx_std_17)
target_link_libraries(foam2umesh
  PUBLIC
  umesh
  PRIVATE
  OpenMP::OpenMP_CXX
  )
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(foam2umesh PRIVATE -Wall -Wno-strict-aliasing)
endif()

option(FOAM2UMESH_STATIC_RUNTIME
  "link libstdc++/libgcc statically into foam2umesh (runs without the GCC module)" OFF)
if (FOAM2UMESH_STATIC_RUNTIME)
  target_link_options(foam2umesh PRIVATE -static-libstdc++ -static-libgcc)
endif()
