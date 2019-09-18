#ifndef __mgmesh12D_h
#define __mgmesh12D_h

// std lib ------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

// conf includes -------------
#include "Printinfo_conf.h"
#include "numeric_vectorM.h"
#include "parallelM.h"
#include "parallel_objectM.h"

// format class
#include <hdf5.h>
// Forwarded class
class MGUtils;
class MGGeomEl;
#define LREF (1)
// ===================================================
/// Class containing information about the mesh
// ===================================================
class MGMesh {  // ==================================================

 public:
  int _dim;                //< Geometrical dimension
  double _mesh2_data[10];  //< Only for coupled mesh

  /// Class utils: directories and parameters
  MGUtils& _mgutils;

  // Multilevel Mesh
  int _NoLevels;  /// Number of Levels (L)

  // processor subdomains (processor P)
  ParallelObjectM _comm;
  int _iproc;     ///< Subdomain (processor P)
  int _n_subdom;  ///< Number of subdomains or processors (subdomain P)

  // Points (only top)
  int* _NoNodes;    ///<  Node for level: index[L+F*_NoLevels]
  double* _xyz;     ///< Coordinates
  double* _xyzo;    ///< Coordinates
  double* _dxdydz;  ///< Coordinates
  int** _node_map;  ///< Node map: connects node level map to top node level ordering: _node_map[L][n]-> top
                    ///< grid
  int* _off_nd[2];  ///< Node map indices subdomain-level offset (_node_map [0]quad [1]lin)

  // Mesh  element topology

  int** _NoElements;   ///< Element level offset (_NoElements[Fem][Level])
  int** _el_map;       ///< Element connectivity (_el_map[Fem][iel + _off_el[Fem][Level+_NoLevels*iproc]])
  int** _el_neighbor;  ///< Element neighbours *(_el_neighbor[Fem][iel + _off_el[Fem][Level+_NoLevels*iproc]])
  int** _off_el;       ///< Element subdomain-level offset (_off_el[Fem][Level+_NoLevels*iproc])
  double* _VolFrac;    ///< Piecewise volume fraction field

  //  Elements
  int _NoFamFEM;        ///< Number of FEM Families (F) Hex, quad, etc
  int* _type_FEM;       ///< Type of VOL+BDRY GeomElement,read from file
  const int _n_GeomEl;  ///< Number of Geometric element classes
  MGGeomEl& _GeomEl;    ///< Geometric element class

  // FEM type  on single element
  int* _ord_FEM;  /// Order (only lagragian)

  // subdomain  Subdomain (no processor P)
  double* _ctrl_dom;  ///< control domain (define a particular subdomain)

  // boundary
  int* _NodeBDgroup;
  double* _normalb;
  double* _dist;  ///< Distance from a wall defined in B_dist
  double* _yplus;

  // ======================================================================================
  // CONSTRUCTOR-DESTRUCTOR
  // ======================================================================================
  // ------------------------------------------------------------------------------------
  ///  Level Mesh constructor
  MGMesh(const ParallelObjectM& comm, MGUtils& mgutils, MGGeomEl& geomel);
  // ------------------------------------------------------------------------------------
  ///  Level Mesh Destructor
  ~MGMesh();

  ///  Substructure destructor
  void clear();
  ///@}
  void idFEMord(const int n);  ///< Set element numbering to identity

  ///@{  \name MESH COORDINATES

  void xcoord(
      double cvect[], const int n_nodes, const int offset) const;  ///< Coordinates (for values x->y->z)

  /// Mesh coordinates (for points (x,y,z))
  void nodesxyz(double xyzvect[], const int n_nodes) const;
  ///@}
  void Translate(const int dir, NumericVectorM& x_old);
  // void MoveMesh(const int dir,NumericVectorM &x_old,int nvars[],int n_fe[]);
  ///@{  \name CONNECTIVITY
  void conn(int gl_conn[], const int ifem, const int indx_mesh) const;  ///< Connectivity for Mesh (indx_mesh)

  /// Connectivity for sub_element (for example hex: hex27->hex20->hex8)
  int sub_conn(int gl_conn[], const int ifem, const int ilev, const int type_el) const;
  ///@}

  // ======================================================================================
  // GET RETURN FUNCTIONS
  // ======================================================================================
  void get_el_nodes(const int el_nnodes, const int bdry, const int Level, const int iel, double xx[])
      const;  ///< Get element nodes
  void get_el_nod_disp(
      const int bdry, const int Level, const int iel,
      /*int  el_conn[],*/ double xx[]) const;
  /// Get element connectivity
  void get_el_conn(const int el_nnodes, const int bdry, const int Level, const int iel, int el_conn[]) const;
  /// Get element neighbours from local processor
  void get_el_neighbor(
      const int el_sides, const int bdry, const int Level, const int iel, int el_neigh[]) const;
  /// Get element neighbours from given processor iproc
  void get_el_neighbor(
      const int el_sides, const int bdry, const int Level, const int iel, int el_neigh[], int iproc) const;
  /// Get element connectivity and node coordinates  from local processor
  void get_el_nod_conn(const int bdry, const int Level, const int iel, int el_conn[], double xx[]) const;
  /// Get element connectivity and node coordinates from given processor iproc
  void get_el_nod_conn(
      const int bdry, const int Level, const int iel, int el_conn[], double xx[], int iproc) const;
  /// Get center element coordinates
  void get_el_ctr(const int el_nnodes, const int bdry, const double* xx_nds, double* el_xm) const;

  void set_node_normal();
  void get_normal_b(
      const double* xx_boundary,  // all surface coordinates <-
      const double x_c[],         // point inside the element <-
      double normal_b[]           // normal+ tang plane ->
      //     double tg0[]    // tang1[0-(DIM-1)]+tang2[DIM-] ->
      ) const;

  void set_node_normal_bc(int bc[], int node_dof[]);
  void set_node_normal_adj_bc(int bc[], int node_dof[]);

  //    void get_adjustment_normal(NumericVectorM &vel_field,int node_dof[],int offset) const ;

  //     void get_normal_b (
  //     const double* xx_boundary,    //< all surface coordinates <-
  //     const double x_c[],           //< point inside the element <-
  //     double* normal_b              //< normal ->
  // ) const ;

  // ======================================================================================
  // PRINT FUNCTIONS
  // ======================================================================================
  // ---------------------------------------------------------------------------
  /// Print volume mesh
  void print(const int Level, const int t_step);
  // ---------------------------------------------------------------------------
  /// Print a subdomain
  void print_subdom_hf5(std::string filename) const;
  // ---------------------------------------------------------------------------
  //   void print_med(std::string filename);           ///< Print mesh in med format
  void print_dist_hf5(std::string filename, double dist[], std::string dir_name) const;
  // ---------------------------------------------------------------------------
  /// Print distance function from the walls
  void print_VolFrac_hf5(std::string filename, std::string dir_name) const;
  // ---------------------------------------------------------------------------
  /// print control domain
  void print_ctrl_dom_hf5(std::string filename, double ctrl_dom[], std::string dir_name) const;
  // ---------------------------------------------------------------------------
  /// Print in xmf format the topology
  void print_xmf_topology(std::ofstream& out, std::string store_file, int nel, int ik);

  /// Print in xmf format the geometry
  void print_xmf_geometry(std::ofstream& out, std::string store_file, int nodes, int ik);
  /// Mesh writing
  void write_c(const int t_step);

  virtual void prova() {}  // to be polymorphic (-> dynamic casting)

  //   void read_hdf5_NPOINTS(std::string   file_name, std::string  setdata_name, int dim_data, int data[]);
  void readf_hdf5_NPOINTS(std::string file_name, std::string setdata_name, int dim_data, int data[]);

 private:
  ///@{  \name FUNCTION FORMAT DEPENDENT -> /CONTRIB/IOFORMAT
  void print_vol(const int Level, const int t_step);  ///< Print volume mesh
  void print_conn_lin_hf5(const int Level);           ///< Print connectivity
  ///@}

  void read_c();  ///< Mesh reading function
};

#endif
