#ifndef __salome_io_h__
#define __salome_io_h__


#define GROUP_NAME_DIR ("/FAS/Mesh_1/ELEME")
#define MESH_NAME_DIR ("/ENS_MAA/Mesh_1/-0000000000000000001-0000000000000000001/MAI/")
#define COORD_NAME_DIR ("/ENS_MAA/Mesh_1/-0000000000000000001-0000000000000000001/NOE/COO")


// Local includes -----------------------------------------------
#include "libmesh_common.h"
#include "mesh_input.h"
#include "mesh_output.h"

// Forward declarations -----------------------------------------
// class MeshBase;
#include "hdf5.h"


namespace libMesh
{
// Forward declarations
 class MeshBase;
// =====================================================================================
// ==================  SALOME.IO class =================================================
// =====================================================================================

// =====================================================================================
// This class implements writing meshes in the .med format.
// For a full description of the med format see doc  *
class SalomeIO : public MeshInput<MeshBase>,public MeshOutput<MeshBase>
{// ====================================================================================
 public:
  // ------------------------------------------------------------------------------------
  // Constructor.  Takes a non-const Mesh reference 
  SalomeIO (
    MeshBase& mesh                  // mesh   <-
  );
  // Constructor.  Takes a reference to a constant mesh object.
  SalomeIO (
    const MeshBase& mesh             // mesh   <-
  );

  // ------------------------------------------------------------------------------------
  // Reads in a mesh in the Dalome *.med format from the ASCII file given by name. 
  virtual void read (
    const std::string& name          // file name   <-
  );
  // ------------------------------------------------------------------------------------
  // This method implements writing a mesh to  the salome *.med format. 
  virtual void write (
    const std::string& /*name*/       // file name   <-
  ){};
  // This method implements writing a mesh with nodal data to a specified file 
  virtual void write_nodal_data (
    const std::string&,
    const std::vector<Number>&,
    const std::vector<std::string>&){};

  // Flag indicating whether or not to write a binary file.
  bool & binary ();

private:
  // Implementation of the read() function.
  virtual void read_mesh (std::istream& in);
  
  // ============================================================================
  
  // ============================================================================
// This function reads the boundary condition and mat groups
void read_bc_mat(
  hid_t file_id,             // file id (hdf5)
  std::string  el_fem_type,  // vol fem salome
  std::string elbd_fem_type,  // boundary fem salome
   int n_elements,
  int n_nodes
//   int * bc_flag,
// int *	      mat_flag,
// std::map<int, int > bdgroup_name // map of boundary names --> id
) ;

  
// This function return the boundary LIBMESH type and the fem Salome names
// (as argument)
int  read_fem_type(  
hid_t file_id,
  std::string & el_fem_type_vol,
  std::string & el_fem_type_bd
);
 
  // This method implements writing a mesh to a specified file.  
  // This will write an ASCII *.msh file.
  virtual void write_mesh (std::ostream& /*out*/){};

  /**
   * This method implements writing a mesh with nodal data to a specified file
   * where the nodal data and variable names are optionally provided.  This
   * will write an ASCII or binary *.pos file, depending on the binary flag.
   */
  void write_post (const std::string&,
                   const std::vector<Number>* = NULL,
                   const std::vector<std::string>* = NULL);
   // void H5Lget_name_by_idx(hid_t arg1, const char* arg2,  arg3,  arg4, int arg5, char[3] arg6, int arg7,  arg8);

  // print
   hid_t print_Ihdf5(hid_t file,const std::string & name, hsize_t dimsf[],int data[]);
  
  /**
   * Flag to write binary data.
   */
  bool _binary;
 


};



// ------------------------------------------------------------
// SalomeIO inline members
inline
SalomeIO::SalomeIO (const MeshBase& mesh) :
  MeshOutput<MeshBase> (mesh),
  _binary        (false)
{
}


inline
SalomeIO::SalomeIO (MeshBase& mesh) :
  MeshInput<MeshBase>  (mesh),
  MeshOutput<MeshBase> (mesh),
  _binary (false)
{}

inline
bool & SalomeIO::binary ()
{
  return _binary;
}

} // namespace libMesh

#endif // #define __salome_io_h__
