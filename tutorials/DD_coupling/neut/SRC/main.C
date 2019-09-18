// libc+++ include ----------------------------------------
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "Solverlib_conf.h"  // Solver options

// include MED library
#ifdef HAVE_MED
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"
#endif
#include <hdf5.h>  // hdf5

#include "FEMUS.h"
#include "InterfaceProjection.h"
#include "MGUtils.h"

// include program
#include "DONDRA.h"  // neutronic interface to donjon-dragon
#include "InterfaceFunctionDD.h"
#include "MeshExtended.h"
#include "mesh_c.h"  // donjon-dragon Cartesian mesh

// ======================================================================================
int main(
    int argc,
    char** argv) {  // ==================================================================================

  argc = argc;
  argv = argv;  // no unused warning

  // MGUtils class (external parameter)
  std::vector<MGUtils*> mgutils;

  // Mesh name files <- mgutils[0] (NUM_MESH=1)
  std::string MeshPosAndName[NUM_MESH];
  for (int i_mesh = 0; i_mesh < NUM_MESH; i_mesh++) {
    mgutils.push_back(new MGUtils(i_mesh + 1));
    MeshPosAndName[i_mesh] = mgutils[i_mesh]->_mesh_dir + mgutils[i_mesh]->_interface_mesh;
    std::cout << " \n P mesh file " << i_mesh + 1 << "= " << MeshPosAndName[i_mesh] << "\n ";
  }
  // main parameter <- mgutils[0]
  int n_steps = stoi(mgutils[0]->_sim_config["nsteps"]);
  double dt = stod(mgutils[0]->_sim_config["dt"]);
  int print_step = stoi(mgutils[0]->_sim_config["printstep"]);
  int itime_0 = stoi(mgutils[0]->_sim_config["itime"]);
  double time = 0.;

  FEMUS P1(*mgutils[0]);  // P1 constructor -> extended system class

  int iproc = P1.get_MGMesh()._iproc;  // set processor

  // setting system -----------------------------------------------------------
  // set system 2 cell vector
  // setting element interfaces  8 -> T  9 ->pwr
  std::vector<int> IDSvec(1);
  IDSvec[0] = 2;
  P1.init_interface(9, IDSvec, 2, mgutils[0]->_interface_mesh, false);  // false -> elem
  P1.init_interface(8, IDSvec, 2, mgutils[0]->_interface_mesh, true);   // true -> point

  // ==========================================================================
  //  system 2 Neutronic Problem DONDRA donjon/dragon
  DONDRA N(iproc);         // neutronic problem
  N.initialize(1., 0.01);  // MW 2700
  N.set_up("IniPowCompo");
  N.set_mesh_med((P1.getMedMesh()));

  int nfa = N.get_n_fuel_array();
  int na = N.get_n_array();
  std::cout << " \n  ==========  Number of fuel array " << nfa << " Number of array " << na << "====== \n";
  // N Problem interfaces. Interface 2=pwr+flux Interface 3=temperature
  N.init_interface(2, 2, MeshPosAndName[0].c_str());
  N.init_interface(3, 2, MeshPosAndName[0].c_str());

  // *****************************************************

  // ==================================================================================
  // time loop neutronics + thermohydraulics
  double bc_value[1];
  N.print_time_xmf(itime_0, n_steps);
  P1.solve_setup(itime_0, time);
  if (iproc == 0) N._mesh_c->print_mesh_med2hdf5("./RESU/mesh_c.h5");

  // ==========================================================================
  //   // From interface 8 Temperature (P1) -> interface 3 (N) at t=0
  //   // Get interface 8 Temperature (P1) -----------------------------------------
  MEDCoupling::MEDCouplingFieldDouble* bdy = NULL;
  bdy = P1.getValuesOnInterface(8, "T", 1);  // take field from P proble
  const MEDCoupling::DataArrayDouble* array_bdy = bdy->getArray();
  // Get interface interface 3 (N) ---------------------------------------------
  InterfaceFunctionDD* fct3 = N.get_interface_fun(3);
  int n_elem = fct3->get_n();
  // setting fuel, coolant  temperature and coolant density field --------------
  MEDCoupling::DataArrayDouble* array = MEDCoupling::DataArrayDouble::New();
  array->alloc(n_elem, 3);
  for (int i_med = 0; i_med < n_elem; i_med++) {
    //     int node_mg   = map_mg[i_med];  // element id in neutronic "mesh"
    //        int node_med   = map_med[i_med];  // element id in neutronic "mesh"
    double tc = array_bdy->getIJ(i_med, 1);
    if (tc < 500) {
      std::cout << i_med << " elemnt " << tc << "T_c bound  500"
                << "\n ";
      tc = 500;
    }
    if (tc > 700) {
      std::cout << i_med << " elemnt " << tc << "T_c bound 700"
                << "\n ";
      tc = 700;
    }
    double tf = tc + 450;
    if (tf > 1150) {
      std::cout << i_med << " elemnt " << tf << "T_f limited 1150"
                << "\n ";
      tf = 1150;
    }
    array->setIJ(i_med, 0, tf);    // Temperature fuel
    array->setIJ(i_med, 1, 0.65);  // Density Coolant
    array->setIJ(i_med, 2, tc);    // Temperature Coolant
  }
  // setting the interface field ----------------------------------------------
  MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
  f->setMesh(fct3->getSupport());  // set mesh
  f->setName("temp");              // set name
  f->setArray(array);              // set array -> f
  f->checkConsistencyLight();      // check f
  N.setBC(3, f);                   // 2-> interface f->field
  N.set_values(3, "NEU", 3, 0);    // 2-> interface 3->n_cmp 0-> first_cmp

  // initial time loop (t=0)
  for (int itime = itime_0; itime < n_steps; itime++) {
    // solve one step problem N
    N.solve_onestep();

    // ==========================================================================
    // From interface 2 pwr (N) -> interface 9 (P1) at t
    //  getting the neutronics output from interface (2-NEU) and set to MEDfield phi
    MEDCoupling::MEDCouplingFieldDouble* phi = NULL;
    phi = N.getValuesPWBUP_elem(2, "NEU");
    phi->setName("power");
    P1.setFieldSource(9, 1, phi);
    // ==========================================================================

    // solve one step problem P1
    P1.solve_onestep(itime_0, itime, print_step, time, dt);

    // print neutronics output
    std::string name_cmp[2];
    name_cmp[0] = "power";
    name_cmp[1] = "burnup";
    std::ostringstream file_name_h5;
    file_name_h5 << "RESU/pw." << itime << ".h5";
    if (iproc == 0) {
      N.print_power_med2hdf5(file_name_h5.str(), *phi, 2, name_cmp);
      N.print_power_xmf(itime, *mgutils[0], file_name_h5.str(), 2, name_cmp);
      std::cout << "Solution at iteration time = " << itime << std::endl;
    }
    // ==========================================================================
    //    From interface 8 Temperature (P1) -> interface 3 (N) at t
    //    Get interface 8 Temperature (P1) -----------------------------------------
    bdy = P1.getValuesOnInterface(8, "T", 1);  // take field from P problem
    const MEDCoupling::DataArrayDouble* array_bdy = bdy->getArray();
    // setting fuel, coolant  temperature and coolant density field --------------
    for (int i_med = 0; i_med < n_elem; i_med++) {
      double tc = array_bdy->getIJ(i_med, 1);
      if (tc < 500) {
        std::cout << i_med << " elemnt " << tc << " T_c bounded 500"
                  << "\n ";
        tc = 500;
      }
      if (tc > 700) {
        std::cout << i_med << " elemnt " << tc << " T_c bounded 700"
                  << "\n ";
        tc = 700;
      }
      double tf = tc + 450;
      if (tf > 1150) {
        std::cout << i_med << " elemnt " << tf << " T_f bounded 1150"
                  << "\n ";
        tf = 1150;
      }
      array->setIJ(i_med, 0, tf);    // Temperature fuel
      array->setIJ(i_med, 1, 0.65);  // Density Coolant
      array->setIJ(i_med, 2, tc);    // Temperature Coolant
    }
    // setting the interface field ----------------------------------------------
    MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
    f->setMesh(fct3->getSupport());  // set mesh
    f->setName("temp");              // set name
    f->setArray(array);              // set array -> f
    f->checkConsistencyLight();      // check f

    N.setBC(3, f);                 // 2-> interface f->field
    N.set_values(3, "NEU", 3, 0);  // 2-> interface 3->n_cmp 0-> first_cmp
    // ==================================================================================

  }  // end itime
  system("rm *.o2m; rm *.l2m");
  array->decrRef();
  return 0;
}
