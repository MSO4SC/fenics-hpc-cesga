// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2006-12-12

#ifndef __SPACE_TIME_FUNCTION_H
#define __SPACE_TIME_FUNCTION_H

#include <dolfin/function/Function.h>

#include <ufc.h>

#include <stdint.h>

namespace dolfin
{
  
  class Mesh;
  
  class SpaceTimeFunction
  {
  public:
    
    /// Create space-time function
    SpaceTimeFunction(Mesh& mesh, Function& Ut);
    
    /// Copy constructor
    SpaceTimeFunction(const SpaceTimeFunction& f);
    
    /// Destructor
    ~SpaceTimeFunction();
    
    /// Evaluate function at time t, giving result in Ut
    void eval(real t);
    
    // Add a space function at time t
    void addPoint(std::string Uname, real t);
    
    // Add a set of functions with arbitrary time steps
    void util_addFiles(std::vector<std::string> filenames);
    
    // Add a set of functions with fixed time step
    void util_addFiles(std::vector<std::string> filenames, real T);
    void util_fileList(std::string basename, int N,
		       std::vector<std::string>& filenames);
    
    /// Return mesh associated with function
    Mesh& mesh();
    
    /// Return interpolant function
    Function& evaluant();
    
  private:
    
    // Pointer to mesh associated with function (null if none)
    Mesh * const mesh_;
    
    // Pointer to evaluant function
    Function * const function_;
    
    // Space functions defining the current time interval (cache)
    Function U0;
    Function U1;
    
    real u0_t;
    real u1_t;
    
    bool u0_t_valid;
    bool u1_t_valid;
    
    std::map<real, std::string> U_files;
    
#ifdef ENABLE_MPIIO
    
    // File headers for DOLFIN's binary file format, repeated here
    // until this information is available outside of DOLFIN
    
    enum Binary_data_t
      {
	BINARY_MESH_DATA,
	BINARY_VECTOR_DATA,
	BINARY_FUNCTION_DATA,
	BINARY_MESH_FUNCTION_DATA
      };
    
    typedef struct
    {
      uint32_t magic;
      uint32_t bendian;
      uint32_t pe_size;
      Binary_data_t type;
    } BinaryFileHeader;
    
    typedef struct
    {
      uint32_t dim;
      uint32_t size;
      real t;
      char name[256];
    } BinaryFunctionHeader;
    
#endif
  };
  
}

#endif
