// #include <dolfin/common/common_includes.h>
// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-11-26
// Last changed: 2006-12-12
//
// Note: this breaks the standard envelope-letter idiom slightly,
// since we call the envelope class from one of the letter classes.

#include <cmath>
#include <iomanip>
#include <sstream>

#include "dolfin/SpaceTimeFunction.h"
#include <dolfin/config/dolfin_config.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/File.h>
#include <dolfin/la/Vector.h>
#include <dolfin/main/MPI.h>

using namespace dolfin;


//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(Mesh& mesh, Function& Ut) :
  mesh_(&mesh),
  function_(&Ut),
  U0(Ut),
  U1(Ut),
  u0_t(0.),
  u1_t(0.),
  u0_t_valid(false),
  u1_t_valid(false)
{
}

//-----------------------------------------------------------------------------
SpaceTimeFunction::SpaceTimeFunction(const SpaceTimeFunction& f) :
  mesh_(f.mesh_),
  function_(f.function_),
  u0_t(f.u0_t),
  u1_t(f.u1_t),
  u0_t_valid(false),
  u1_t_valid(false)
{
}

//-----------------------------------------------------------------------------
SpaceTimeFunction::~SpaceTimeFunction()
{

}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::eval(real t)
{
  
  std::map<real, std::string>::iterator it1;
  std::map<real, std::string>::iterator it0;
  
  if (U_files.size() == 0)
  {
    error("Cannot interpolate on zero sample files");
  }
  
  // NOTE: t is the current time in the primal referential t \in [0,primal_Tend]
  // Find element in U_files so that element < t
  
  // Select it1 such that the time t1 is just after t
  it1 = U_files.upper_bound(t);
  
  // If t == T, we need to step back one
  if (it1 == U_files.end())
  {
    --it1;
  }
  
  it0 = it1;
  --it0;
  
  real t0 = (*it0).first;
  real t1 = (*it1).first;
  
  if (t0 != t0 || t1 != t1)
  {
    error("At least one of the iteration times"
	  "used for interpolation is a Nan.");
  }
  
  std::string name0 = (*it0).second;
  std::string name1 = (*it1).second;
  
  if (t0 != u0_t || !u0_t_valid)
  {
    File file0(name0);
    u0_t_valid = true;
    u0_t = t0;
    file0 >> U0.vector();
  }
  
  if (t1 != u1_t || !u1_t_valid)
  {
    File file1(name1);
    u1_t_valid = true;
    u1_t = t1;
    file1 >> U1.vector();
  }
  
  // Compute weights (linear Lagrange interpolation)
  real w0 = (t1 - t) / (t1 - t0);
  real w1 = (t - t0) / (t1 - t0);
  
  cout << "S0: t = " << t0 << "; name0 = " << name0 << "; w0 = " << w0 << endl;
  cout << "S1: t = " << t1 << "; name1 = " << name1 << "; w1 = " << w1 << endl;
  
  // Compute interpolated value
  evaluant().vector() = 0.0;
  evaluant().vector().axpy(w0, U0.vector());
  evaluant().vector().axpy(w1, U1.vector());
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::addPoint(std::string Uname, real t)
{
  U_files[t] = Uname;
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::util_addFiles(std::vector<std::string> filenames)
{
  
#ifdef ENABLE_MPIIO
  
  int counter = 0;
  int num_files = filenames.size();
  
  for (std::vector<std::string>::iterator it = filenames.begin();
       it != filenames.end(); ++it)
  {
    std::string filename = *it;
    
    MPI_File fh;
    MPI_Offset byte_offset;
    BinaryFileHeader hdr;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE,
		      MPI_STATUS_IGNORE);
    
    byte_offset = sizeof(BinaryFileHeader);

    uint nfunc;
    MPI_File_read_at_all(fh, byte_offset, &nfunc, sizeof(uint), MPI_BYTE,
			 MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);
    BinaryFunctionHeader f_hdr;
    MPI_File_read_at_all(fh, byte_offset, &f_hdr, 
			 sizeof(BinaryFunctionHeader),
			 MPI_BYTE, MPI_STATUS_IGNORE);
    
    // Temporary load function, and parse time stamp
    addPoint(filename, f_hdr.t);
    
    MPI_File_close(&fh);
    
    counter++;
  }
#else
  error("MPI I/O required for space time functions with arbitrary time step");
#endif
  
}

//-----------------------------------------------------------------------------
void SpaceTimeFunction::util_addFiles(std::vector<std::string> filenames,
                                      real T)
{
  //FIXME: For now we assume a fixed time step
  
  int counter = 0;
  int num_files = filenames.size();
  
  if (num_files == 1)
  {
    error("Divide by zero");
  }
  for (std::vector<std::string>::iterator it = filenames.begin();
       it != filenames.end(); ++it)
  {
    std::string filename = *it;
    
    // OK guys this is *only* valid if we do the right thing i.e:
    // - num_files is the number of samples
    // - T is the measure of the time interval for solving the dual problem
    //	 i.e [sampling_start_time, primal_end_time]
    real t = T * real(counter) / real(num_files - 1);
    std::cout << "add intermediate time t = " << t << std::endl;
    
    addPoint(filename, t);
    
    counter++;
  }
  if (counter == 0)
  {
    error("Counter irremediably stayed stuck at zero.");
  }
  for (std::map<real, std::string>::const_iterator it = U_files.begin();
       it != U_files.end(); ++it)
  {
    std::cout << std::setw(4) << it->first << " : " << it->second << std::endl;
  }
  
}
//-----------------------------------------------------------------------------
void SpaceTimeFunction::util_fileList(std::string basename, int N,
                                      std::vector<std::string>& filenames)
{
  filenames.clear();
  
  // OK so N is the number of samples spanning [T0,T1]
  // Then there are N - 1 time intervals
  for (int sample_id = 0; sample_id < N; ++sample_id)
  {
    std::stringstream filename, number;
    number.fill('0');
    number.width(6);
    
    number << sample_id;
    
    filename << basename;
    filename << number.str();
#ifdef ENABLE_MPIIO
    filename << ".bin";
#else
    filename << "_" << MPI::processNumber() << ".bin";
#endif
    filename << std::ends;
    
    filenames.push_back(filename.str());
  }
  if (filenames.size() == 0)
  {
    error("Trying to interpolate over zero samples");
  }
}

//-----------------------------------------------------------------------------
Mesh& SpaceTimeFunction::mesh()
{
  dolfin_assert(mesh_);
  return *mesh_;
}
//-----------------------------------------------------------------------------
Function& SpaceTimeFunction::evaluant()
{
  dolfin_assert(function_);
  return *function_;
}
//-----------------------------------------------------------------------------



