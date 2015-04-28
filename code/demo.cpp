/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2012, Marianna Madry
*  All rights reserved.
*
*  Contact: marianna.madry@gmail.com
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of Willow Garage, Inc. nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
********************************************************************/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "libsvm_file.h"
#include "string_kernel.h"

using std::cerr;
using std::endl;
using std::cout;
using std::string;
using std::vector;


void usage(const char *exec_name) {
  cerr << "---------------------------------------------\n"
          "Compute a string kernel and save it to a file\n"
          "in the libSVM kernel format\n"
          "By Marianna Madry (marianna.madry@gmail.com)\n"
          "---------------------------------------------\n";
  cerr << "Usage: " << exec_name << " <kernel_file>\n";
  cerr << "Args:\n"
          "  <kernel_file>: output file to save kernel values\n"
          "                 in libsvm format\n";
}

int main(int argc, char **argv) {
  if (argc != 2) {
    usage(argv[0]);
    exit(1);
  }
  string kernel_file(argv[1]);

  // Kernel parameters
  const float c = 1e12;
  const int normalize = 1;
  const int symbol_size = 255;  // A size of an alphabet
  const int max_length = 1000;  // A maximum sequence length
  int kn = 2;                   // A level of susbsequence matching
  double lambda = 0.5;          // A decay factor

  // Prepare dummy data
  vector<string> dummy_data;
  dummy_data.push_back("The idea behind digital computers may be explained by saying that these machines are intended to carry out any operations which could be done by a human. Alan Turing."); // An original quote from Alan Turing
  dummy_data.push_back("The idea may be explained by saying that these intelligent machines are intended to carry out all tasks which could be done by a human."); // A changed sentence
  dummy_data.push_back("AAGCTAGCTAGCAAGCTAGCTAGC"); // An example of a DNA sequence 
  dummy_data.push_back("TAGTAGCTAAAGCTAGCTTTA");

  // Prepare labels for dummy data
  vector<int> dummy_labels;
  dummy_labels.push_back(1);
  dummy_labels.push_back(1);
  dummy_labels.push_back(-1);
  dummy_labels.push_back(-1);

  // Main computations
  StringKernel<float> string_kernel(c, normalize, symbol_size, max_length, kn, lambda);
  string_kernel.set_data(dummy_data);
  string_kernel.compute_kernel();

  // Save kernel to file
  if (write_libsvm_kernel(kernel_file, dummy_labels, string_kernel))
    std::cout << "Kernel saved in the libsvm format to: " << kernel_file << std::endl;
  else
  {
    std::cout << "Error: Cannot write to file: " << kernel_file << std::endl;
    exit(1);
  }
}
