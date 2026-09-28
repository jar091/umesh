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
# BHAC (MPI-AMRVAC GRMHD) snapshot on a spherical MKS grid -> N
# spatially coherent .umesh parts, plus diagnostics (see README.md)
# ------------------------------------------------------------------
add_executable(bhac2umesh
  bhac2umesh.cpp
  BhacSnapshot.cpp
  )
target_compile_features(bhac2umesh PRIVATE cxx_std_17)
target_link_libraries(bhac2umesh
  PUBLIC
  umesh
  )

# reloads .umesh parts and checks indices/NaNs/orientation/volumes
add_executable(bhacUmeshCheck
  bhacUmeshCheck.cpp
  )
target_compile_features(bhacUmeshCheck PRIVATE cxx_std_17)
target_link_libraries(bhacUmeshCheck
  PUBLIC
  umesh
  )

# CPU emulation of barney's unstructured-volume path
add_executable(barneyVolumeEmu
  barneyVolumeEmu.cpp
  )
target_compile_features(barneyVolumeEmu PRIVATE cxx_std_17)
target_link_libraries(barneyVolumeEmu
  PUBLIC
  umesh
  )

# OpenMP only speeds up the diagnostic modes (--ref-render, barneyVolumeEmu)
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
  target_link_libraries(bhac2umesh PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(barneyVolumeEmu PRIVATE OpenMP::OpenMP_CXX)
endif()
