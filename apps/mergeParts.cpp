// ======================================================================== //
// Copyright 2018-2026 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

/* merges several umesh files (e.g., neighbouring parts of a data-parallel
   decomposition) into one, so that a mesh split into N parts can be
   rendered on N/k GPUs. Vertices shared between the inputs are _not_
   merged (see mergeMeshes()); elements and per-vertex scalars are
   appended with their vertex indices shifted. */

#include "umesh/UMesh.h"

namespace umesh {

  void usage(const std::string error="")
  {
    if (error != "")
      std::cerr << "\nError : " << error  << "\n\n";

    std::cout << "Usage: ./umeshMergeParts -o <out.umesh> <in0.umesh> <in1.umesh> ...\n\n";
    exit(error != "");
  };

  extern "C" int main(int ac, char **av)
  {
    std::string outFileName;
    std::vector<std::string> inFileNames;
    for (int i=1;i<ac;i++) {
      const std::string arg = av[i];
      if (arg == "-h")
        usage();
      else if (arg == "-o")
        outFileName = av[++i];
      else if (arg[0] != '-')
        inFileNames.push_back(arg);
      else
        usage("unknown cmd-line arg '"+arg+"'");
    }

    if (outFileName == "") usage("no output file specified");
    if (inFileNames.empty()) usage("no input files specified");

    std::vector<UMesh::SP> inputs;
    for (auto fileName : inFileNames) {
      std::cout << "loading umesh from " << fileName << std::endl;
      inputs.push_back(UMesh::loadFrom(fileName));
      if (!inputs.back()->perVertex)
        throw std::runtime_error(fileName+" has no per-vertex scalars");
    }
    UMesh::SP merged = mergeMeshes(inputs);
    inputs.clear();
    std::cout << "merged " << inFileNames.size() << " parts: "
              << merged->toString() << std::endl;
    std::cout << "writing to " << outFileName << std::endl;
    merged->saveTo(outFileName);
  }

} // ::umesh
