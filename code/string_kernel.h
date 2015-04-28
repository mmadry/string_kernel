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

#ifndef _STRING_KERNEL_H_
#define _STRING_KERNEL_H_

#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include "data_set.h"

template<class k_type>
class StringKernel {
 public:
  /** Constructor, sets kernel parameters. */
  StringKernel(const float c, const int normalize, const int symbol_size,
               const size_t max_length, int kn, double lambda)
      : _c(c), _normalize(normalize), _symbol_size(symbol_size),
        _max_length(max_length), _kn(kn), _lambda(lambda),
        _string_data(0), _kernel(0)
        {}

  ~StringKernel()
  {
    for (size_t i = 0; i < _string_data->size(); i++)
      delete[] _kernel[i];
    delete [] _kernel;
    delete _string_data;
  }

  /** Set the dataset to be used by the kernel. */
  void set_data(const std::vector<std::string> &strings);

  /** Calculate the kernel. */
  void compute_kernel();

  /** Return pointer to kernel matrix. */
  k_type **values() const
  {
    assert(_kernel);
    return _kernel;
  }

  /** Return the size of the NxN kernel. */
  size_t size() const
  {
    assert(_string_data);
    return _string_data->size();
  }

 protected:
  const float _c;
  const int _normalize;
  const int _symbol_size;
  const size_t _max_length;
  const int _kn;
  const double _lambda;
  DataSet *_string_data;
  k_type **_kernel;

 private:
  k_type kernel(const DataElement &x, const DataElement &y) const;
  void run_kernel_dp(const std::vector<k_type> &norms, k_type **K) const;
};


template<class k_type>
void StringKernel<k_type>::set_data(const std::vector<std::string> &strings)
{
  assert(strings.size() > 0);
  _string_data = new DataSet(_max_length, _symbol_size);
  _string_data->load_strings(strings);
}


template<class k_type>
k_type StringKernel<k_type>::kernel(const DataElement &x, const DataElement &y) const
{
  k_type **Kd[2];
  for (int i = 0; i < 2; i++)
  {
    Kd[i] = new k_type *[x.length + 1];
    for (int j = 0; j < x.length + 1; j++)
    {
      Kd[i][j] = new k_type[y.length + 1];
    }
  }

  // Dynamic programming
  // Initialize the two two-dimensional Kd-arrays
  // [one for the current i (i%2), one for i-1 ((i+1)%2) ]
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < (x.length + 1); j++)
    {
      for (int k = 0; k < (y.length + 1); k++)
      {
        Kd[i][j][k] = (i + 1) % 2;
        // we start with i==1, i.e. K'_i-1(.,.)==1.
        //        If we increase i to 2, we switch
      }
    }
  }

  // Calculate all the Kd and Kdd
  for (int i = 1; i <= (_kn - 1); i++)
  {
    // Set the Kd's to nought for those lengths of s and t
    // where s (or t) has exactly length i-1 and t (or s)
    // has length >=i-1  this part of the matrix has the shape
    // of a capital letter L upside down
    for (int j = (i - 1); j <= (x.length - 1); j++)
    {
      Kd[i % 2][j][i - 1] = 0;
    }
    for (int j = (i - 1); j <= (y.length - 1); j++)
    {
      Kd[i % 2][i - 1][j] = 0;
    }
    for (int j = i; j <= (x.length - 1); j++)
    {
      k_type Kdd = 0;
      for (int m = i; m <= (y.length - 1); m++)
      {
        if (x.attributes[j - 1] != y.attributes[m - 1])
        {
          // ((.))-1 is because indices start with 0 (not with 1)
          Kdd = _lambda * Kdd;
        }
        else
        {
          Kdd = _lambda * (Kdd + (_lambda * Kd[(i + 1) % 2][j - 1][m - 1]));
        }
        Kd[i % 2][j][m] = _lambda * Kd[i % 2][j - 1][m] + Kdd;
      }
    }
  }

  // Calculation of K
  k_type sum = 0;
  for (int i = _kn; i <= x.length; i++)
  {
    for (int j = _kn; j <= y.length; j++)
    {
      if (x.attributes[((i)) - 1] == y.attributes[((j)) - 1])
      {
        // ((.))-1 is because indices start with 0 (not with 1)
        sum += _lambda * _lambda * Kd[(_kn - 1) % 2][i - 1][j - 1];
      }
    }
  }

  // Delete
  for (int j = 0; j < 2; j++)
  {
    for (int i = 0; i < x.length + 1; i++)
    {
      delete[] Kd[j][i];
    }
  }
  for (int i = 0; i < 2; i++)
  {
    delete[] Kd[i];
  }

  return sum;
}

template <class k_type>
void StringKernel<k_type>::run_kernel_dp(const std::vector<k_type> &norms,
                                         k_type **K) const
{
  assert(_string_data);

  for (size_t i = 0; i < _string_data->size(); i++)
  {
    for (size_t j = 0; j < _string_data->size(); j++)
    {
      if (K[j][i] == -1)
      {
        K[j][i] = kernel(_string_data->elements()[j], _string_data->elements()[i]);

        if (_normalize)
          K[j][i] /= sqrt(norms[i] * norms[j]);

        K[i][j] = K[j][i];
      }
    }
  }
}


template<class k_type>
void StringKernel<k_type>::compute_kernel()
{
  assert(_string_data);

  // Initialize kernel, array of floats
  _kernel = new k_type *[_string_data->size()];
  for (size_t i = 0; i < _string_data->size(); i++)
    _kernel[i] = new k_type[_string_data->size()];

  // start with all K filled with -1's, then only calculate kernels as needed
  for (size_t i = 0; i < _string_data->size(); i++)
    for (size_t j = 0; j < _string_data->size(); j++)
      _kernel[i][j] = -1;


  // get values for normalization, it is computed for elements in diagonal
  std::vector<k_type> norms(_string_data->size());
  if (_normalize)
  {
    for (size_t i = 0; i < _string_data->size(); i++)
    {
      norms[i] = kernel(_string_data->elements()[i],
                        _string_data->elements()[i]);
      _kernel[i][i] = 1;
    }
  }

  // Compute kernel using dynamic programming
  run_kernel_dp(norms, _kernel);
}



#endif
