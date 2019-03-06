#ifndef __meshc_HXX__
#define __meshc_HXX__

#include "Solverlib_conf.h"


#include <string>



  class MeshC {
  // DATA -------------------------------------------------------------------------------
  private:
  int _nx; int _ny; int _nz;
  int _hx; int _hy; int _hz;
  public:
  
 
  double * _x;
  double * _y;
  double * _z;
    
  // Constructor-Destructor--------------------------------------------------------------
  MeshC(){};
  MeshC(
    int nex_in,int ney_in,int nez_in // elements
  ){
    _nx=nex_in;_ny=ney_in;_nz=nez_in;
     _hx=0.;_hy=0.;_hz=0.;   // Cartesian dimension domain Vol=hx hy hz
     int n_nodes_mesh_c=(_nx+1)*(_ny+1)*(_nz+1);
    _x=new double[n_nodes_mesh_c];
    _y=new double[n_nodes_mesh_c];
    _z=new double[n_nodes_mesh_c];
  };
   MeshC(int nex_in,int ney_in,int nez_in,
	 int hx_in,int hy_in,int hz_in  // Cartesian dimension domain Vol=hx hy hz
	){
    _nx=nex_in;_ny=ney_in;_nz=nez_in;
     _hx=hx_in;_hy=hy_in;_hz=hz_in;  // Cartesian dimension domain Vol=hx hy hz
     int n_nodes_mesh_c=(_nx+1)*(_ny+1)*(_nz+1);
    _x=new double[n_nodes_mesh_c];
    _y=new double[n_nodes_mesh_c];
    _z=new double[n_nodes_mesh_c];
  };
  ~MeshC(){return;}
  
  double & get_x(){return *_x;}
  double & get_y(){return *_y;}
  double & get_z(){return *_y;}
  int get_nex() const {return _nx;}
  int get_ney() const {return _ny;}
  int get_nez() const {return _nz;}
  int get_hx() const {return _hx;}
  int get_hy() const {return _hy;}
  int get_hz() const {return _hz;}
  int get_n_nodes(){return (_nx+1)*(_ny+1)*(_nz+1);}
  int get_n_elements(){return _nx*_ny*_nz;}
  void  print_mesh_med2hdf5(
  std::string file_name                  ///< file name (where to print)
  );
  
  };

#endif
