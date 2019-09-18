#include "DONDRA.h"
#include "mesh_c.h"
using namespace boost;
using namespace std;
using namespace ganlib;

#include <hdf5.h>

#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionDD.h"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"
#endif
#include <iomanip>
#include "MGUtils.h"
// using namespace MEDCoupling;

//=======================================================================================
/// This funciton is a constructor fro the class DONDRA
DONDRA::DONDRA() : _mesh_c(NULL) {
  // =========================================================================
  cout << "New DONDRA object constructed.'" << endl;
}
//=======================================================================================
/// This funciton is a constructor fro the class DONDRA
DONDRA::DONDRA(int iproc) : _iproc(iproc), _mesh_c(NULL) {
  // =========================================================================
  cout << "New DONDRA object constructed.'" << endl;
}
//=======================================================================================
/// This function is the initialization function for the class  DONDRA
void DONDRA::initialize(
    double power,   ///< thermal reactor power in MW
    double Lref     ///<  reference Length
) {                 // =================================================================================
  _power_ = power;  // thermal reactor power in MW
  _Lref = Lref;
}
//=======================================================================================
/// This function is the set up function for the data exchange for the class  DONDRA
void DONDRA::set_up(std::string filename_in  ///< name init power  (ex IniPowCompo_.c2m)
) {
  // construct the Lifo stack for IniPowCompo -------------------------------------------
  cout << "DONDRA::set_up" << endl;
  LifoPtr ipLifo1 = LifoPtr(new Lifo());
  ipLifo1->pushEmpty("Fmap", "LCM");
  ipLifo1->pushEmpty("Matex", "LCM");
  ipLifo1->pushEmpty("Cpo", "LCM");
  //   ipLifo1->pushEmpty("CpRefl", "LCM");
  ipLifo1->pushEmpty("Track", "LCM");

  //   std::ostringstream filename;
  // call IniPowCompo Cle-2000 procedure ------------------------------------------------
  //   filename<< filename_in;// << _iproc;
  //   Cle2000Ptr IniPowCompo = Cle2000Ptr(new Cle2000(filename.str().c_str(), 0, ipLifo1));

  Cle2000Ptr IniPowCompo = Cle2000Ptr(new Cle2000("IniPowCompo", 2, ipLifo1));
  IniPowCompo->exec();
  cout << "IniPowCompo execution completed" << endl;

  // recover the output LCM objects -----------------------------------------------------
  ipLifo1->node("Fmap", _Fmap);
  ipLifo1->node("Matex", _Matex);
  ipLifo1->node("Cpo", _Cpo);
  //   ipLifo1->node("CpRefl", _CpRefl);
  ipLifo1->node("Track", _Track);
  IntPtrConst stateVector = _Fmap->getInt("STATE-VECTOR");
  _n_fuel_array = stateVector[0] * stateVector[1];
  long npar = stateVector[7];

  // empty the Lifo stack
  while (ipLifo1->getMax() > 0) ipLifo1->pop();

  // Getting the mesh from Cartesian mesh (_mesh_c) -------------------------------------
  IntPtrConst stateVector_matex0 = _Matex->getInt("STATE-VECTOR");

  const int nx = stateVector_matex0[7]; /* _mesh_c._nx=nx;*/
  const int ny = stateVector_matex0[8]; /*_mesh_c._ny=ny;*/
  const int nz = stateVector_matex0[9]; /* _mesh_c._nz=nz;*/
  _mesh_c = new MeshC(nx, ny, nz);
  _n_array = (nx) * (ny) * (nz);
  int n_nodes_mesh_c = (nx + 1) * (ny + 1) * (nz + 1);
  //   _mesh_c._xxx = new double [n_nodes_mesh_c] ;
  //   _mesh_c._yyy = new double [n_nodes_mesh_c];
  //   _mesh_c._zzz = new double [n_nodes_mesh_c];
  FloatPtrConst meshxTab = _Matex->getFloat("MESHX");
  FloatPtrConst meshyTab = _Matex->getFloat("MESHY");
  FloatPtrConst meshzTab = _Matex->getFloat("MESHZ");
  for (int k = 0; k <= nz; k++) {
    for (int j = 0; j <= ny; j++) {
      for (int i = 0; i <= nx; i++) {
        const int ind = k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
        _mesh_c->_x[ind] = meshxTab[i] * _Lref;
        _mesh_c->_y[ind] = meshyTab[j] * _Lref;
        _mesh_c->_z[ind] = meshzTab[k] * _Lref;
      }
    }
  }
  // Open an existing file and dataset --------------------------------------------------
  hid_t file = H5Fcreate("../RESU/coord.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t dimsf[2];
  dimsf[0] = n_nodes_mesh_c;
  dimsf[1] = 1;
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);
  hid_t dtsetx = H5Dcreate(file, "X1", H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dtsety = H5Dcreate(file, "X2", H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t dtsetz = H5Dcreate(file, "X3", H5T_NATIVE_DOUBLE, dtsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the data to the dataset using default memory space, file space, and transfer properties.
  hid_t status = H5Dwrite(dtsetx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_mesh_c->_x[0]);
  status = H5Dwrite(dtsety, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_mesh_c->_y[0]);
  status = H5Dwrite(dtsetz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_mesh_c->_z[0]);

  // clean ------------------------------------------------------------------------------
  H5Dclose(dtsetx);
  H5Dclose(dtsety);
  H5Dclose(dtsetz);
  H5Fclose(file);
  ipLifo2 = LifoPtr(new Lifo());
  //   PowComponent = Cle2000Ptr(new Cle2000("PowComponent", 0, ipLifo2));

  std::ostringstream filename_pow;
  // call IniPowCompo Cle-2000 procedure ------------------------------------------------
  filename_pow << "PowComponent";  // << _iproc;
  PowComponent = Cle2000Ptr(new Cle2000(filename_pow.str().c_str(), 0, ipLifo2));
  return;
}

// ======================================================================================
void DONDRA::solve_onestep() {
  // iteration loop
  ClcmPtr Flux, Thm;
  float_32 densB = 2000.;
  float Keff_conv = 1.0;
  int iter = 1;
  //    int continueLoop = 1;

  //    ipLifo2 = LifoPtr(new Lifo());
  //     PowComponent = Cle2000Ptr(new Cle2000("PowComponent", 0, ipLifo2));
  // ipLifo2 = LifoPtr(new Lifo());
  //   PowComponent = Cle2000Ptr(new Cle2000("PowComponent", 0, ipLifo2));

  std::ostringstream filename_pow;
  // call IniPowCompo Cle-2000 procedure ------------------------------------------------
  //   filename_pow<<"PowComponent";// << _iproc;
  //   PowComponent = Cle2000Ptr(new Cle2000(filename_pow.str().c_str(), 0, ipLifo2));

  const int nx = _mesh_c->get_nex();
  const int ny = _mesh_c->get_ney();
  const int nz = _mesh_c->get_nez();

  cout << "DONDRA: ITERATION NUMBER:" << iter << " iteration = " << iter << endl;

  // construct the Lifo stack for PowComponent
  ipLifo2->push("Fmap", _Fmap);
  ipLifo2->push("Matex", _Matex);

  if (iter == 1) {
    ipLifo2->pushEmpty("Flux", "LCM"); /*ipLifo2->pushEmpty("MacroP", "LCM");*/
  } else {
    ipLifo2->push("Flux", Flux); /*ipLifo2->push("MacroP", _MacroP);*/
  }
  ipLifo2->push("Cpo", _Cpo);
  ipLifo2->push("Track", _Track);
  ipLifo2->push("iter", iter);
  ipLifo2->push("powi", float_32(_power_));
  ipLifo2->push("densB", densB);

  // call PowComponent Cle-2000 procedure
  cout << "call PowComponent->exec" << endl;
  PowComponent->exec();
  cout << "PowComponent execution completed" << endl;
  ipLifo2->node("Flux", Flux);
  FloatPtrConst Keff = Flux->getFloat("K-EFFECTIVE");
  cout << "DONDRA: iter=" << iter << " ------------- Keffective=" << Keff[0] << endl;
  Keff_conv = Keff[0];

  //   IntPtrConst stateVector_matex = Matex->getInt("STATE-VECTOR");
  //   cout << "\n STATE-VECTOR Matex " << endl;
  //   for(int i=0; i<40; i++) {
  //     cout <<  stateVector_matex[i] << " ";
  //   }
  //
  //   cout << "\n powerTab with length= " << _mylength << endl;
  //   for(int i=0; i<_mylength; i++) {
  //     cout <<  powerTab[i] << " ";
  //   }
  // //      communicator_.send(iter, "powerTab", mylength, powerTab);
  //
  //
  //
  //   IntPtrConst stateVector_br = Fmap->getInt("STATE-VECTOR");
  //   cout << "\n STATE-VECTOR Fmap new" << endl;
  //   for(int i=0; i<40; i++) {
  //     cout <<  stateVector_br[i] << " ";
  //   }
  //   FloatPtrConst irradiationTab = Fmap->getFloat("BURN-INST");
  //   cout << "\n irradiationTab with length= " <<_mylength << endl;
  //   for(int i=0; i<_mylength; i++) {
  //     cout <<  powerTab[i] << "";
  //   }
  //      communicator_.send(iter, "irradiationTab", mylength, irradiationTab);

  //      // receive thermo-hydraulics information
  //      ClcmPtr Jpmap = Fmap->getClcm("PARAM");
  //      IntPtr myIntPtr(new int_32[1]); myIntPtr[0] = 2;
  //      for (int ipar=0; ipar<npar; ipar++) {
  //        ClcmPtr Kpmap = Jpmap->getClcm(ipar);
  //        StringPtrConst pname = Kpmap->getString("P-NAME");
  //        if (*pname == "T-FUEL      ") {
  //          FloatPtr myArray(new float_32[mylength]);
  // //          communicator_.recv(iter, "fuel_temperature", mylength, myArray);
  //          Kpmap->put("P-VALUE", myArray, mylength);
  //          Kpmap->put("P-TYPE", myIntPtr, 1);
  //        } else if (*pname == "D-COOL      ") {
  //          FloatPtr myArray(new float_32[mylength]);
  // //          communicator_.recv(iter, "water_density", mylength, myArray);
  //          Kpmap->put("P-VALUE", myArray, mylength);
  //          Kpmap->put("P-TYPE", myIntPtr, 1);
  //        } else if (*pname == "T-COOL      ") {
  //          FloatPtr myArray(new float_32[mylength]);
  // //          communicator_.recv(iter, "water_temperature", mylength, myArray);
  //          Kpmap->put("P-VALUE", myArray, mylength);
  //          Kpmap->put("P-TYPE", myIntPtr, 1);
  //        }
  //      }
  _Fmap->val();

  // empty the Lifo stack
  while (ipLifo2->getMax() > 0) ipLifo2->pop();

  // receive the convergence flag
  //      communicator_.recv(iter, "continueLoop", continueLoop);
  cout << "DONDRA: Value of continueLoop : " << iter << " at iteration " << iter << endl;
  //    }
  //   cout << "DONDRA: close the Calcium communicator" << endl ;
  // //    communicator_.terminate();
  cout.precision(10);
  cout << "DONDRA: converged K-effective=" << Keff_conv << endl;

  return;
}

// ========================================================================
// This function gets the  mesh name from the med file (meshNames)
// and sets the id "interface_id" and the mesh to the interfaces vector
// (_interfaceFunMap) through the index "interface_name" of the
// interface-functions. The interface is  an element map only.
void DONDRA::init_interface(
    const int interface_name, int interface_id,
    const std::string& medfile_name  // medfile name    (in)
) {  // =========================================================================

  ostringstream name_id;
  name_id << interface_id;
  std::vector<std::string> vG(1);
  std::string id_name = name_id.str().c_str();
  vG[0] = name_id.str().c_str();

  // Reading mesh Names from med file  ----------------------------------------
  //   std::string mesh_dir=_mg_utils->_mesh_dir;
  //   std::string localFile=mesh_dir+medfile_name;

  std::vector<std::string> meshNames = MEDCoupling::GetMeshNames(medfile_name.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if (meshNames.size() < 1) { std::cout << " FEMUS::setMesh : no meshes in the file'"; }
  //   std::string localMeshName = meshNames[0];//index_medmesh=0
  //   // Reading group names
  //   std::vector<std::string> GroupNames =
  //     MEDLoader::GetMeshGroupsNames(medfile_name.c_str(), localMeshName.c_str());
  //   int nGroups0 = GroupNames.size();

  // extra=================================================
  // get the med_mesh (support) -------------------------------------------------
  MEDCoupling::MEDCouplingUMesh* support;
  int id_level = -1;
  if (interface_id < 10) { id_level = 0; }
  support = MEDCoupling::ReadUMeshFromGroups(medfile_name.c_str(), meshNames[0].c_str(), id_level, vG);
  support->zipCoords();
  std::cout << "FEMUS::setInterfaces: read mesh from interface " << interface_id << " with name "
            << interface_name << " at level " << id_level << "\n";

  int nElement_med = support->getNumberOfCells();                // n MED-elements (_n now is element)
  const MEDCoupling::DataArrayDouble* d = support->getCoords();  // Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples();                      //  MED-nodes
  assert(n_nodes_med == support->getNumberOfNodes());            // check  MED-nodes
  int dim_med = d->getNumberOfComponents();                      //  MED dimension
  int n_nod_el = (dim_med == 3) ? 8 : 4;                         // for Cartesian cell only
  double* xyz_med = new double[nElement_med * dim_med];          // MED el centers

  // Computing the MED center (xyz_med) -------------------------------------
  {
    std::vector<int> nodes1;  // element nodes
    double x_m[3];            // tmp buffer

    for (int ielem = 0; ielem < nElement_med; ielem++) {
      for (int idim = 0; idim < dim_med; idim++) x_m[idim] = 0.;  // zeros
      support->getNodeIdsOfCell(ielem, nodes1);                   // element nodes
      for (int inode = 0; inode < n_nod_el; inode++) {
        for (int im = 0; im < dim_med; im++) { x_m[im] += d->getIJ(nodes1[inode], im); }

      }                // end inode
      nodes1.clear();  // clear element node vector

      for (int im = 0; im < dim_med; im++) {
        xyz_med[dim_med * ielem + im] = x_m[im] / n_nod_el;
        std::cout << ielem << " x" << im << "= " << x_m[im] / n_nod_el << " ";
      }
      std::cout << "\n";
    }
  }
  //===========================================================

  // Create the interface
  InterfaceFunctionDD* fun = new InterfaceFunctionDD;
  // set the mesh interface
  fun->set_mesh_interface_elemID(*_mesh_c, support, interface_name);
  // setting the fun in the map _interfaceFunMap at bd_name_id
  _interfaceFunMap[interface_name] = fun;  // added to the map

  return;
}

// ======================================================================================
// This function gets the  the value of the variable with id number "variable_id" on
// nodes on the boundary with identity "id":   variable_name -> "BUND-PW"
MEDCoupling::MEDCouplingFieldDouble* DONDRA::getValuesPW_elem(
    int id,                    // int boundary identity   (in)
    const char* variable_name  // system name             (in)
) {  // ==================================================================================
  // from LibMesh function fct ----------------------------------------------------------
  InterfaceFunctionDD* fct = get_interface_fun(id);
  if (fct == NULL) return NULL;
  int nElem = fct->getSupport()->getNumberOfCells();
  int n_elem = fct->get_n();

  int* map_mg = fct->get_map_mg();  // map_mg[i_med]: Med-Mesh-> DonJon-Mesh
  int* map_med = fct->get_map_med();
  const int nx = _mesh_c->get_nex();
  const int ny = _mesh_c->get_ney();
  const int nz = _mesh_c->get_nez();
  int n_element_mesh_c = (nx + 1) * (ny + 1) * (nz + 1);

  // ------------------------------------------------------------------------------------
  // power from the neutronic system -> send reactor physics information
  IntPtrConst indexfmTab = _Fmap->getInt("BMIX");  // index assemblies(indexfmTab !=0)
  FloatPtrConst powerTab = _Fmap->getFloat("BUND-PW");
  double* pwr = new double[n_element_mesh_c];
  int counter_ipower = 0;
  //   double val[2];
  for (int k = 0; k < n_element_mesh_c; k++) {
    double val = 0;
    if (indexfmTab[k] != 0) {
      val = powerTab[counter_ipower];
      counter_ipower++;
    }
    pwr[k] = val;
  }

  // array data vector (array) to fill the field f --------------------------------------
  MEDCoupling::DataArrayDouble* array = MEDCoupling::DataArrayDouble::New();
  array->alloc(n_elem, 1);
  for (int i_med = 0; i_med < n_elem; i_med++) {  // filling array
    int node_mg = map_mg[i_med];                  // element id in neutronic "mesh"
    array->setIJ(node_mg, 0, pwr[node_mg]);
  }

  // field f (ON_CELLS) of the interface  fct -------------------------------------------
  MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
  f->setMesh(fct->getSupport());  // set Med-mesh
  f->setName(variable_name);      // set name
  f->setArray(array);             // set array -> f

  // check and clean --------------------------------------------------------------------
  f->checkConsistencyLight();  // check f
                               //   fct->set_field(f);
                               //   fct->printOn(std::cout,id);   // print fct
  array->decrRef();            // delete array
  delete[] pwr;

  // Print -------------------------------------------------------------------------------
  std::cout << "\n getValuesPW_elem with " << counter_ipower << " fuel elements over " << nElem
            << "total elements" << std::endl;

  return f;
}

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
// ------------------------------
//   variable_name ="BUND-PW"   -> power
//   variable_name ="BURN-INST" -> burn-up
// -------------------------------
MEDCoupling::MEDCouplingFieldDouble* DONDRA::getValuesPWBUP_elem(
    int id,                    // int boundary identity   (in)
    const char* variable_name  // system name             (in)
) {                            // ========================================================================
  // from LibMesh function fct
  InterfaceFunctionDD* fct = get_interface_fun(id);
  if (fct == NULL) return NULL;
  int nElem_med = fct->getSupport()->getNumberOfCells();
  int n_elem_med = fct->get_n();
  int n_element_mesh_c = _mesh_c->get_n_elements();
  int* map_mg = fct->get_map_mg();
  int* map_med = fct->get_map_med();

  // ------------------------------------------------------------------------------------
  // take power from the neutronic system and send reactor physics information
  IntPtrConst indexfmTab = _Fmap->getInt("BMIX");               // get assembly reactor matrix
  FloatPtrConst powerTab = _Fmap->getFloat("BUND-PW");          // get power
  FloatPtrConst irradiationTab = _Fmap->getFloat("BURN-INST");  // get  burn-up

  // data into pwr (power) and irr (burn-up) vector
  double* pwr = new double[n_element_mesh_c];  // power
  double* irr = new double[n_element_mesh_c];  // burn-up
  double val[2];
  int counter_ipower = 0;
  double potmax = 0.;
  for (int k = 0; k < n_element_mesh_c; k++) {
    if (indexfmTab[k] != 0) {
      pwr[k] = powerTab[counter_ipower];
      potmax += pwr[k];
      irr[k] = irradiationTab[counter_ipower];
      counter_ipower++;
    } else {
      pwr[k] = -1.e-4;
      irr[k] = -1.e-4;
    }
  }
  double potmax1 = 0.;
  for (int k = 0; k < n_element_mesh_c; k++) {
    if (indexfmTab[k] != 0) {
      pwr[k] = pwr[k] * counter_ipower / potmax;
      potmax1 += pwr[k];
    }
  }

  std::cout << " ======================= " << potmax1 << "  potmax1 \n";
  // array function to fill f -----------------------------------------------------------
  MEDCoupling::DataArrayDouble* array = MEDCoupling::DataArrayDouble::New();
  array->alloc(n_elem_med, 2);
  // filling array
  for (int i_med = 0; i_med < n_elem_med; i_med++) {
    int node_mg = map_mg[i_med];           // element id in neutronic "mesh"
    array->setIJ(i_med, 0, pwr[node_mg]);  // power
    array->setIJ(i_med, 1, irr[node_mg]);  // burn-up
  }

  // field f ----------------------------------------------------------------------------
  MEDCoupling::MEDCouplingFieldDouble* f = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS);
  f->setMesh(fct->getSupport());  // set mesh
  f->setName("power");            // set name
  f->setArray(array);             // set array -> f

  // check and  clean -------------------------------------------------------------------
  f->checkConsistencyLight();  // check f
  //   fct->set_field(f); fct->printOn(std::cout,id);   // print fct
  array->decrRef();  // delete array
  delete[] pwr;
  delete[] irr;

  // print ------------------------------------------------------------------------------
  std::cout << "\n GetValuesPWBup_elem: " << counter_ipower << " fuel elements (" << nElem_med << " total)"
            << std::endl;
  return f;
}

// =========================================================================
void DONDRA::set_values(
    int id_boundary_name,       ///< identity interface name (in)
    std::string mgsystem_name,  ///< system name             (in)
    int n_cmp,                  ///<  variable system    (in)
    int first_cmp               ///< from variable system      (in)

) {  // ======================================================================

  InterfaceFunctionDD* fct = _interfaceFunMap[id_boundary_name];  // interface-function
  if (fct == NULL) return;

  const MEDCoupling::MEDCouplingUMesh* support = fct->getSupport();
  int nNodes_med = support->getNumberOfNodes();
  int n_nodes_mg = fct->get_n();
  //   int * map_mg  = fct->get_map_mg();
  //   int * map_med = fct->get_map_med();
  // Cartesian grid for neutronics
  int ne_x = _mesh_c->get_nex();
  int ne_y = _mesh_c->get_ney();
  int ne_z = _mesh_c->get_nez();

  int nElement_med = support->getNumberOfCells();                // n MED-elements (_n now is element)
  const MEDCoupling::DataArrayDouble* d = support->getCoords();  // Med-mesh coordinates
  int n_nodes_med = d->getNumberOfTuples();                      //  MED-nodes
  assert(n_nodes_med == support->getNumberOfNodes());            // check  MED-nodes
  int dim_med = d->getNumberOfComponents();                      //  MED dimension
  int n_nod_el = (dim_med == 3) ? 8 : 4;                         // for Cartesian cell only
  double* xyz_med = new double[nElement_med * dim_med];          // MED el centers

  // Computing the MED center (xyz_med) -------------------------------------
  {
    std::vector<int> nodes1;  // element nodes
                              //     double x_m[3];//tmp buffer

    for (int ielem = 0; ielem < nElement_med; ielem++) {
      //       for(int idim=0; idim<dim_med; idim++) x_m[idim]=0.; // zeros
      support->getNodeIdsOfCell(ielem, nodes1);  // element nodes
      for (int im = 0; im < dim_med; im++) {
        double sum = 0.;
        for (int inode = 0; inode < n_nod_el; inode++) { sum += d->getIJ(nodes1[inode], im); }
        xyz_med[dim_med * ielem + im] = sum / n_nod_el;
      }                // end inode
      nodes1.clear();  // clear element node vector
    }
  }

  double bc_value[3];
  IntPtrConst stateVector = _Fmap->getInt("STATE-VECTOR");
  IntPtrConst indexfmTab = _Fmap->getInt("BMIX");  // get assembly reactor matrix
                                                   //   _n_fuel_array = stateVector[0] * stateVector[1];
  long npar = stateVector[7];

  int n_element_mesh_c = _mesh_c->get_n_elements();

  //  map mesh_c -> fuel
  int* map11 = new int[n_element_mesh_c];
  int counter_ipower = 0;
  for (int k = 0; k < n_element_mesh_c; k++) {
    if (indexfmTab[k] != 0) {
      map11[k] = counter_ipower;
      counter_ipower++;
    } else
      map11[k] = -1;
  }
  assert(counter_ipower == _n_fuel_array);

  // receive thermo-hydraulics information
  FloatPtr myArray_tf(new float_32[_n_fuel_array]);
  FloatPtr myArray_td(new float_32[_n_fuel_array]);
  FloatPtr myArray_tc(new float_32[_n_fuel_array]);
  for (int i_mg = 0; i_mg < nElement_med; i_mg++) {
    int index_2 = (int)(xyz_med[3 * i_mg + 0] * ne_x / 3.66) +
                  ((int)(xyz_med[3 * i_mg + 1] * ne_y / 3.66)) * ne_x +
                  ((int)(xyz_med[3 * i_mg + 2] * ne_z / 6.67)) * ne_y * ne_x;
    //     int gl_node_bd = map_mg[i_mg];  // mg  node
    if (indexfmTab[index_2] != 0) {
      //     int node_med  = map_med[i_mg];  // med node
      fct->eval(i_mg, 3, bc_value);
      //
      int indx54 = map11[index_2];
      if (indx54 < 0 || indx54 > _n_fuel_array) {
        std::cout << " ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR ";
      }
      myArray_tf[indx54] = bc_value[0];
      myArray_td[indx54] = bc_value[1];
      myArray_tc[indx54] = bc_value[2];
    }
  }

  ClcmPtr Jpmap = _Fmap->getClcm("PARAM");
  IntPtr myIntPtr(new int_32[1]);
  myIntPtr[0] = 2;
  for (int ipar = 0; ipar < npar; ipar++) {
    ClcmPtr Kpmap = Jpmap->getClcm(ipar);
    StringPtrConst pname = Kpmap->getString("P-NAME");
    if (*pname == "C-BORE      ") continue;
    IntPtrConst ptype = Kpmap->getInt("P-TYPE");
    //     if (ptype[0] != 2) throw Cle2000Exception("THM failure");
    if (*pname == "T-FUEL      ") {  // fuel_temperature
      Kpmap->put("P-VALUE", myArray_tf, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    } else if (*pname == "D-COOL      ") {  //  water_density
      Kpmap->put("P-VALUE", myArray_td, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    } else if (*pname == "T-COOL      ") {  // water_temperature
      Kpmap->put("P-VALUE", myArray_tc, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    }
  }

  //       _Fmap->val();

  // empty the Lifo stack
  //    while(ipLifo2->getMax() > 0) ipLifo2->pop();

  //     for(int jdim=first_cmp; jdim<first_cmp+n_cmp; jdim++) {
  //       const int kdof_top = mgsyst->_node_dof[Level][gl_node_bd+jdim*offset];
  //       old_sol_top->set(kdof_top , bc_value[jdim-first_cmp]);
  //     }

  delete[] xyz_med;
  //   delete []map_mg;
  //   delete []map_med;
  delete[] map11;
  return;
}

// =========================================================================
void DONDRA::set_tc_values(
    const MEDCoupling::MEDCouplingFieldDouble& bdy,
    int id_boundary_name,       ///< identity interface name (in)
    std::string mgsystem_name,  ///< system name             (in)
    int itime) {                // ======================================================================

  InterfaceFunctionDD* fct = _interfaceFunMap[id_boundary_name];  // interface-function
  if (fct == NULL) return;

  const MEDCoupling::MEDCouplingUMesh* support = fct->getSupport();
  //   int nNodes_med = support->getNumberOfNodes();
  int n_nodes_mg = fct->get_n();
  int* map_mg = fct->get_map_mg();
  int* map_med = fct->get_map_med();

  const MEDCoupling::DataArrayDouble* array_bdy = bdy.getArray();
  //   MGSolBase * mgsyst=get_eqs(mgsystem_name.c_str());
  //   int Level=_mg_mesh->_NoLevels-1;
  //   int offset= _mg_mesh->_NoNodes[Level];
  //   NumericVectorM *old_sol_top= mgsyst->x_old[Level];
  //
  //
  //   MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
  //   std::vector<int> mat=ext_mesh->_mat_id;
  //   std::vector<int> bc_id=ext_mesh->_bc_id;
  //
  //   std::vector<double> src_value(3);
  //   double *myArray=new double[n_nodes_mg];
  double bc_value[3];
  IntPtrConst stateVector = _Fmap->getInt("STATE-VECTOR");
  //   _n_fuel_array = stateVector[0] * stateVector[1];
  long npar = stateVector[7];

  // receive thermo-hydraulics information
  ClcmPtr Jpmap = _Fmap->getClcm("PARAM");
  IntPtr myIntPtr(new int_32[1]);
  myIntPtr[0] = 2;

  FloatPtr myArray_tf(new float_32[_n_fuel_array]);
  FloatPtr myArray_td(new float_32[_n_fuel_array]);
  FloatPtr myArray_tc(new float_32[_n_fuel_array]);
  double max = 0.;
  double min = 100000000.;
  for (int i_mg = 0; i_mg < _n_fuel_array; i_mg++) {
    //      int gl_node_bd = map_mg[i_mg];  // mg  node
    int node_med = map_med[i_mg];  // med node
    double tc = array_bdy->getIJ(node_med, 2);
    double tf = array_bdy->getIJ(node_med, 0);
    double densc = array_bdy->getIJ(node_med, 1);
    //     tc *=1.10;
    //        if(tc>max) max=tc;if(min>tc) min=tc;
    //     if(tc< 500)  {
    //       std::cout << i_mg << " elemnt "<< tc<< " coolant T_c increased to 500" <<"\n ";
    //       tc= 500;
    //     }
    //     if(tc> 700)  {
    //       std::cout << i_mg << " elemnt "<<tc << " coolant T_c element limited 700" <<"\n ";
    //       tc= 700;
    //     }
    //     double tf=(tc+300)*1.5;
    //     if(tf > 1150) {
    //       std::cout << i_mg << " elemnt "<<  tf<< " fuel element T_f limited to  1150" <<"\n ";
    //       tf=1150;
    //
    //     }
    //    fct->eval(node_med,3, bc_value);
    std::cout << " ****** " << node_med << " " << tc << " " << tf << "\n";

    myArray_tf[node_med] = tf;
    myArray_td[node_med] = 0.65;
    myArray_tc[node_med] = tc;
  }

  for (int ipar = 0; ipar < npar; ipar++) {
    ClcmPtr Kpmap = Jpmap->getClcm(ipar);
    StringPtrConst pname = Kpmap->getString("P-NAME");
    if (*pname == "C-BORE      ") continue;
    IntPtrConst ptype = Kpmap->getInt("P-TYPE");
    //     if (ptype[0] != 2) throw Cle2000Exception("THM failure");

    if (*pname == "T-FUEL      ") {  // fuel_temperature
      //        std::cout << "\n\n\n --------------- min max T ------" << max << " ------- " << min << "\n\n
      //        \n";
      Kpmap->put("P-VALUE", myArray_tf, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    } else if (*pname == "D-COOL      ") {  //  water_density
      Kpmap->put("P-VALUE", myArray_td, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    } else if (*pname == "T-COOL      ") {  // water_temperature
      Kpmap->put("P-VALUE", myArray_tc, _n_fuel_array);
      Kpmap->put("P-TYPE", myIntPtr, 1);
    }
  }

  //       _Fmap->val();

  // empty the Lifo stack
  //   while(ipLifo2->getMax() > 0) ipLifo2->pop();

  //     for(int jdim=first_cmp; jdim<first_cmp+n_cmp; jdim++) {
  //       const int kdof_top = mgsyst->_node_dof[Level][gl_node_bd+jdim*offset];
  //       old_sol_top->set(kdof_top , bc_value[jdim-first_cmp]);
  //     }

  return;
}

// =======================================================================
void DONDRA::setBC(
    int name,  ///< interface name
               //   int  n_cmp,            ///< number of components
    const MEDCoupling::MEDCouplingFieldDouble*
        field) {  // =======================================================================

  InterfaceFunctionDD* fct = get_interface_fun(name);  // interface-function
  if (fct == NULL) return;
  fct->set_field(field);
  // #ifdef PRINT_MED
  //   fct->printOn(std::cout,name);
  // #endif

  return;
}
// ======================================================================================
/// This function prints in hdf format on linear cartesian grids (HEX27->HEX8)
/// the following power quantites: n_cmp=1 -> power; n_cmp=2 power+burnup
void DONDRA::print_power_med2hdf5(
    std::string file_name,                     ///< file name (where to print)
    MEDCoupling::MEDCouplingFieldDouble& phi,  ///< mesh field
    int n_cmp,                                 ///< stampa n_cmp=1->power;n_cmp=2 power+burnup
    std::string name_cmp[]                     ///< dataset name vector  ("power","burnup")
) {  // ===================================================================================

  // setup ------------------------------------------------------------------------------
  int n_element_mesh_c = _mesh_c->get_n_elements();
  // take power from the neutronic system and send reactor physics information
  IntPtrConst indexfmTab = _Fmap->getInt("BMIX");               // get assembly reactor matrix
  FloatPtrConst powerTab = _Fmap->getFloat("BUND-PW");          // get power
  FloatPtrConst fluxTab = _Fmap->getFloat("FLUX-AV");           // get power
  FloatPtrConst irradiationTab = _Fmap->getFloat("BURN-INST");  // get  burn-up

  // data into pwr (power) and irr (burn-up) vector
  double* pwr = new double[n_element_mesh_c];   // power
  double* irr = new double[n_element_mesh_c];   // burn-up
  double* flux = new double[n_element_mesh_c];  // burn-up
  double val[2];
  int counter_ipower = 0;
  double potmax = 0.;
  for (int k = 0; k < n_element_mesh_c; k++) {
    if (indexfmTab[k] != 0) {
      pwr[k] = powerTab[counter_ipower];
      potmax += pwr[k];
      irr[k] = irradiationTab[counter_ipower];
      flux[k] = fluxTab[counter_ipower];
      counter_ipower++;
    } else {
      pwr[k] = -1.e-12;
      irr[k] = -1.e-12;
      flux[k] = 1.e-12;
    }
  }
  double potmax1 = 0.;
  for (int k = 0; k < n_element_mesh_c; k++) {
    if (indexfmTab[k] != 0) {
      pwr[k] = pwr[k] * counter_ipower / potmax;
      potmax1 += pwr[k];
    }
  }

  // printing in the hdf5 file  ---------------------------------------------------------
  hid_t file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  // file
  hsize_t dimsf[2];
  dimsf[0] = n_element_mesh_c;
  dimsf[1] = 1;  // dataspace dimension vector dimsf

  // power
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace dtsp
  hid_t dtsetxp = H5Dcreate(
      file, "power", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  hid_t status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pwr[0]);
  H5Dclose(dtsetxp);  // close dataspace
  // burn up
  dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace dtsp
  dtsetxp = H5Dcreate(
      file, "burnup", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &irr[0]);
  H5Dclose(dtsetxp);  // close dataspace

  // flux
  dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace dtsp
  dtsetxp = H5Dcreate(
      file, "flux", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &flux[0]);
  H5Dclose(dtsetxp);  // close dataspace

  // clean an close ---------------------------------------------------------------------
  delete[] pwr;
  delete[] irr;
  H5Fclose(file);

  return;
}

// ======================================================================================
void DONDRA::print_power_xmf(
    const int t_step, MGUtils& mgutils, std::string file_name_h5,
    int n_cmp,              ///< stampa n_cmp=1->power;n_cmp=2 power+burnup
    std::string name_cmp[]  ///< dataset name vector  ("power","burnup")
) {
  //  setup -----------------------------------------------------------------------------
  const int ndigits = 5;  //_mgutils.get_par("ndigits");
  std::string elem_type("Hexahedron");
  int n_node_elem = 8;
  // 3D setup
  int nxyz[3];
  nxyz[0] = _mesh_c->get_nex();
  nxyz[1] = _mesh_c->get_ney();
  nxyz[2] = _mesh_c->get_nez();
  // 2D setup
  if (DIMENSION == 2) {
    elem_type = "Quadrangular";
    n_node_elem = 4;
    nxyz[2] = 0;
  }
  //  Mesh ------------------------------------------------------------------------------
  //   const MGMesh& mgmesh=_mgmesh;
  int n_nodes_lin = (nxyz[0] + 1) * (nxyz[1] + 1) * (nxyz[2] + 1);  //
  std::cout << "DONDRA::print_power_xmf -> nx ny nz " << nxyz[0] << " " << nxyz[1] << " " << nxyz[2]
            << std::endl;
  int n_elements_lin = _mesh_c->get_n_elements();  //
  std::cout << "DONDRA::print_power_xmf -> n_nodes " << n_nodes_lin << "  n_elements  " << n_elements_lin
            << std::endl;
  //
  //   // get parameters
  std::string inout_dir = mgutils._inout_dir;
  std::string basesol = "donjsol";
  std::string basemesh = mgutils.get_file("BASEMESH");
  std::string contrib_dir = mgutils.get_file("CONTRIB_DIR");
  std::string aux_xdmf = mgutils.get_file("AUX_XDMF");
  std::string connlin = mgutils.get_file("CONNLIN");
  //   const double dt = mgutils.get_par("dt");
  double dt = stod(mgutils._sim_config["dt"]);
  //
  //   // files
  // //   std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
  // //   conn_file /* <<femus_dir  <<  "/" << appl_dir << "/"  << myapp << "/" << input_dir */ << basemesh;
  std::ostringstream topol_file;  // topology file file (mesh_conn_lin.h5)
  topol_file << basemesh << connlin << ".h5";
  // //   conn_file << ".h5";
  std::ostringstream coord_time_file;  // connectivity file (mesh_conn_lin.h5)
  coord_time_file /*<<femus_dir  << "/"<< output_dir << "/"*/ << basemesh << "." << std::setw(ndigits)
                                                              << std::setfill('0') << t_step << ".h5";

  std::ostringstream filename;  //  solution file xmf
  filename << inout_dir << basesol << "." << setw(ndigits) << setfill('0') << t_step;
  std::ostringstream attr_file;  //  solution file h5
  /* attr_file << filename.str()<< ".h5"; */
  attr_file << basesol << "." << setw(ndigits) << setfill('0') << t_step << ".h5";
  filename << ".xmf";

  // solution file  xmf
  std::ofstream out(filename.str().c_str());
  std::cout << " " << filename.str().c_str() << std::endl;
  // //  ++++++++++++ Header ++++++++++++++
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<Xdmf> \n"
      << "<Domain> \n";
  out << "<Grid Name=\"Mesh\"> \n";  // =======================================================================
  // time
  out << "<Time Value =\"" << t_step * dt << "\" /> \n";
  // +++++ Topology ++++++++++
  out << " <Topology TopologyType= \"" << elem_type.c_str() << "\"  Dimensions=\"" << n_elements_lin
      << "\"> \n";
  out << " <DataItem Dimensions= \"" << n_elements_lin << " " << n_node_elem
      << "\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\"> \n";
  out << "mesh_c.h5:MESH0CONN \n";
  out << " </DataItem> \n";
  out << " </Topology> \n";
  // +++++++  Geometry +++++++++++++++++
  out << " <Geometry GeometryType=\"X_Y_Z\"> \n";
  out << " <DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_nodes_lin
      << " 1\" Format=\"HDF\"> ";
  out << "mesh_c.h5:X1 \n";
  out << " </DataStructure> \n";
  out << " <DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_nodes_lin
      << " 1\" Format=\"HDF\"> ";
  out << "mesh_c.h5:X2 \n";
  out << " </DataStructure> \n";
  out << " <DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_nodes_lin
      << " 1\" Format=\"HDF\"> ";
  out << "mesh_c.h5:X3 \n";
  out << " </DataStructure> \n";
  out << " </Geometry> \n";
  //   // ++++  Attributes ++++++++++++
  out << "<Attribute Name=\"pw\" AttributeType=\"Scalar\" Center=\"Cell\"> \n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements_lin
      << " 1\" Format=\"HDF\">  \n";
  out << "pw." << t_step << ".h5:power \n";
  out << "</DataItem>\n";
  out << "</Attribute>\n";
  out << "<Attribute Name=\"bup\" AttributeType=\"Scalar\" Center=\"Cell\"> \n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements_lin
      << " 1\" Format=\"HDF\">  \n";
  out << "pw." << t_step << ".h5:burnup \n";
  out << "</DataItem>\n";
  out << "</Attribute>\n";
  out << "<Attribute Name=\"flux\" AttributeType=\"Scalar\" Center=\"Cell\"> \n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements_lin
      << " 1\" Format=\"HDF\">  \n";
  out << "pw." << t_step << ".h5:flux \n";
  out << "</DataItem>\n";
  out << "</Attribute>\n";
  out << "</Grid>\n";  // =================================================================
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();

  return;
}

// =====================================
void DONDRA::print_time_xmf(const int t_in, const int t_end) {
  const int ndigits = 5;
  std::ostringstream filename;  //  solution file xmf
  filename << "./RESU/donjsol.time.xmf";
  std::ofstream out(filename.str().c_str());

  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"/homesd/msandro/software/femus//contrib/Xdmf.dtd\"[]> \n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  out << "<Grid Name=\"sol.msh1\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
  for (int it = t_in; it < t_end; it++) {
    out << "<xi:include href=\"donjsol." << setw(ndigits) << setfill('0') << it << ".xmf"
        << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid[1])\" > \n";
    out << "<xi:fallback /> \n";
    out << " </xi:include> \n";
  }
  out << "</Grid> \n";
  out << "</Domain> \n";
  out << "</Xdmf> \n";

  return;
}

void DONDRA::print_med(std::string namefile, const MEDCoupling::MEDCouplingFieldDouble& field /* Level1*/) {
  //  int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
  //   int n_elements = _mgmesh._NoElements[0][_NoLevels-1];
  //   int nodes_el = _mgmesh._type_FEM[0];
  //
  //   double * coord; coord = new double[n_nodes*3];
  //
  //   for(int i=0; i<n_nodes; i++) {
  //     coord[i*3+1]=0.; coord[i*3+2]=0.;
  //     for(int idim=0; idim<_mgmesh._dim; idim++) {
  //       coord[i*3+idim]=_mgmesh._xyz[i+idim*n_nodes];
  // //        std::cout << " "<< coord[i*3+idim];
  //     }
  // //      std::cout<< " "<< coord[i*3+2] << std::endl;
  //   }
  //
  //   std::cout << " " << n_nodes << " " << n_elements << " " << nodes_el << std::endl;
  //   int  Level=_mgmesh._NoLevels-1;
  //   int icount =0;
  //
  //   int * conn; conn=new int [n_elements*nodes_el];
  //   for(int  iproc = 0; iproc <_mgmesh._n_subdom; iproc++) {
  //     for(int el = _mgmesh._off_el[0][iproc*_mgmesh._NoLevels+Level];
  //          el <_mgmesh._off_el[0][iproc*_mgmesh._NoLevels+Level+1]; el++) {
  //         for(int  i = 0; i < nodes_el; i++) {
  //          conn[icount] = _mgmesh._el_map[0][el*nodes_el+i]; icount++;
  //       }
  //     }
  //   }
  //
  // _med_mesh
  //  MEDCoupling::MEDCouplingUMesh
  //   MEDCoupling::MEDCouplingUMesh *mesh=MEDCoupling::MEDCouplingUMesh::New("Mesh_1",_mgmesh._dim);
  //   mesh->allocateCells(n_elements);
  //   for(int  i = 0; i < n_elements; i++) mesh->insertNextCell(MED_EL_TYPE,nodes_el,conn+i*nodes_el);
  //   mesh->finishInsertingCells();
  //
  //   MEDCoupling::DataArrayDouble *coordarr=MEDCoupling::DataArrayDouble::New();
  //   coordarr->alloc(n_nodes,3);
  //   std::copy(coord,coord+n_nodes*3,coordarr->getPointer());
  //   mesh->setCoords(coordarr);
  //   MEDLoader::WriteUMesh(namefile.c_str(), mesh, true);
  // //   coordarr->decRef();
  // //   mesh->decRef();
  //
  //   delete[] coord;

  //   MEDLoader::WriteUMesh(namefile.c_str(), _med_mesh, true);

  //  set up
  //   const int offset=_mgmesh._NoNodes[_NoLevels-1];
  //   double* sol=new double[n_nodes+1];
  //
  //   MEDCoupling::MEDCouplingFieldDouble * f  =
  //   MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES); f->setMesh(mesh);
  //   f->setName(_var_names[0].c_str()); MEDCoupling::DataArrayDouble *array =
  //   MEDCoupling::DataArrayDouble::New(); array -> alloc(n_nodes, 1);
  //
  //   // print quad -------------------------------------
  //   for(int ivar=0; ivar<1; ivar++)        {
  //     std::string var_name = _var_names[ivar];
  //     for(int i=0; i< n_nodes; i++) {
  //       sol[i]  = (*x_old[Level])(_node_dof[Level][i+ivar*offset])*_refvalue[ivar];
  //       array ->setIJ(i,0,sol[i]);
  //     }
  //   }
  //
  //   f->setArray(array);

  std::string s = "./RESU/my_test3_sol_1";
  s += ".med";
  MEDCoupling::WriteField(namefile.c_str(), &field, false);
  return;
}

// ===========================================================================================
// ===========================================================================================

// ======================================================================================
