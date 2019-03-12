// include local class
#include "MGUtils.h"
#include "EquationsMap.h"
#include "FEMUS.h"
#include "InterfaceProjection.h"
#include "IbUtils.h"

#include "Domain_conf.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"

// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

MEDCoupling::MEDCouplingFieldDouble  *GetFieldComponent( MEDCoupling::MEDCouplingFieldDouble* fieldToExtract, int compID);

// =======================================
// MAIN PROGRAM FOR NAVIER STOKES TESTS
// =======================================

int main ( int argc, char ** argv ) {

// //     // MED FIELDS
  MEDCoupling::MEDCouplingFieldDouble * PW_VolumeFraction, *NodeVolumeFraction;
//
  std::vector<MGUtils *> mgutils;
  std::string MeshPosAndName[NUM_MESH];

  for ( int i_mesh = 0; i_mesh < NUM_MESH; i_mesh++ ) {
      mgutils.push_back ( new MGUtils ( i_mesh + 1 ) );
      MeshPosAndName[i_mesh] = mgutils[i_mesh]->_mesh_dir + mgutils[i_mesh]->_interface_mesh;
      std::cout << " \n P mesh file " << i_mesh + 1 << "= " << MeshPosAndName[i_mesh] << "\n ";
      }

  int    n_steps = stoi ( mgutils[0]->_sim_config["nsteps"] );
  double      dt = stod ( mgutils[0]->_sim_config["dt"] );
  int print_step = stoi ( mgutils[0]->_sim_config["printstep"] );
  int    itime_0 = stoi ( mgutils[0]->_sim_config["itime"] );
  int    restart = stoi ( mgutils[0]->_sim_config["restart"] );
  int    Levels = (int) ( mgutils[0]->_geometry["nolevels"] );
  double time    = 0.;

  // CONSTRUCTION OF FEMUS PROBLEM ------------------------------------------
  FEMUS P ( *mgutils[0] ); //  parameter list <- mgutils[0]
  P.solve_setup ( itime_0, time );                // initial time loop (t=0)
  
  mgutils[0]->read_temp("DATA/Equations.in");
  
  const int NavierStokesSolver = stoi(mgutils[0]->_temp_info["MG_NavierStokes"]);
  const int Proc = P.get_proc();
  std::vector<int> Group;
  Group.push_back ( 1 );
  P.init_interface ( 11, Group, 2, mgutils[0]->_interface_mesh, true );

  const MEDCoupling::MEDCouplingUMesh * FemusMesh = P.getUMesh ( 11 );

  MEDCoupling::MEDCouplingUMesh * PalaMesh = MEDCoupling::ReadUMeshFromFile ( "MESH/body2.med" );
  MEDCoupling::MEDCouplingFieldDouble * FemVel = NULL;
  
if(NavierStokesSolver==1){
  FemVel = P.getValuesOnInterface ( 11, "NS0", DIMENSION, 0 );
}
else{

  const MEDCoupling::MEDCouplingFieldDouble * FemU = P.getValuesOnInterface ( 11, "NS0X", 1, 0 );
  const MEDCoupling::MEDCouplingFieldDouble * FemV = P.getValuesOnInterface ( 11, "NS0Y", 1, 0 );

  FemVel = MEDCoupling::MEDCouplingFieldDouble::MeldFields ( FemU, FemV );
#if DIMENSION==3
  FemVel = MEDCoupling::MEDCouplingFieldDouble::MeldFields ( FemVel, FemW );
  const MEDCoupling::MEDCouplingFieldDouble * FemW = P.getValuesOnInterface ( 11, "NS0Z", 1, 0 );
#endif
}
  
  IbUtils IBmover = IbUtils ( PalaMesh, FemusMesh, Proc );
  std::vector<double> Vel ( DIMENSION );
  Vel[0] = 0;
  Vel[1] = 0.;
  Vel[DIMENSION-1] = 0.;
  
  IBmover.InitSolidVelocity ( Vel );
  double Translation[3] = {0, 0, 0};
  IBmover.InitPosition ( Translation, FemVel );
  FemVel = IBmover.GetInterpolatedSolidVelocity();

  P.setFieldSource ( 11, 2, FemVel );
  
if(NavierStokesSolver==1){
      P.setFieldSource ( 11, 2, FemVel );
      P.write_Boundary_value ( 11, "NS0", DIMENSION, 0 );
}
else
{
      P.setFieldSource ( 11, 2, GetFieldComponent(FemVel,0) );
      P.write_Boundary_value ( 11, "NS0X", 1, 0 );
      P.setFieldSource ( 11, 2, GetFieldComponent(FemVel,1) );
      P.write_Boundary_value ( 11, "NS0Y", 1, 0 );
#if DIMENSION==3
      P.setFieldSource ( 11, 2, GetFieldComponent(FemVel,2) );
      P.write_Boundary_value ( 11, "NS0Z", 1, 0 );
#endif
}

  FemVel->decrRef();
  
  
  NodeVolumeFraction = IBmover.GetFluidNodeColor();
  P.setFieldSource ( 11, 1, NodeVolumeFraction );
  P.write_Boundary_value ( 11, "IB2", 1, 0 );
  P.write_Boundary_value ( 11, "IB1", 1, 0 );

  mgutils[0]->set_IB_info ( &IBmover );

  PW_VolumeFraction = IBmover.GetFluidCellColor();
  P.SetPieceFieldOnYdist ( PW_VolumeFraction );
  P.setMeshName ( mgutils[0]->_interface_mesh );

  // INITIALIZATION OF EQUATIONS TO SOLVE -----------------------------------

  P.solve_onestep ( itime_0, 0, print_step, time, dt, IB_VOL, IB_VOL ); // solving P

  // INTERFACES FOR PERIODIC CHANNEL REALIZATION
  std::vector<int> Group_in;
  Group_in.push_back ( 16 );
  P.init_interface ( 16, Group_in, 2, mgutils[0]->_interface_mesh, true );
  
  std::vector<int> Group_out;
  Group_out.push_back ( 14 );
  P.init_interface ( 14, Group_out, 2, mgutils[0]->_interface_mesh, true );
  
  MEDCoupling::MEDCouplingUMesh * OutletMesh = P.getUMesh ( 14 )->deepCopy();
  MEDCoupling::MEDCouplingUMesh * InletMesh = P.getUMesh ( 16 )->deepCopy();
  double DeltaOverlap[3] = {-0.41, 0, 0};
  OutletMesh->translate (DeltaOverlap);
  
  InterfaceProjection *PeriodicProjector = new InterfaceProjection();
  PeriodicProjector->FillParameters(OutletMesh, InletMesh, Boundary);
  
  for ( int itime = itime_0 + 1; itime <= itime_0 + n_steps; itime++ ) {
      time += dt;
      P.solve_onestep ( itime_0, itime, print_step, time, dt, 0, 19 ); // solving P.
      
      MEDCoupling::MEDCouplingFieldDouble *OutVel = P.getValuesOnInterface ( 14, "NS0X", 1, 0 );
      MEDCoupling::MEDCouplingFieldDouble *InVel  = PeriodicProjector->InterpolatedField(OutVel);
      P.setFieldSource(16, 2, InVel);
      P.write_Boundary_value(16, "NS0X", 1, 0);
      OutVel->decrRef();
      InVel->decrRef();
      }

      
  PW_VolumeFraction->decrRef();   
  NodeVolumeFraction->decrRef(); 
  PalaMesh->decrRef(); 
  FemusMesh->decrRef(); 
  // CLEANING MEMORY --------------------------------------------------------
  P.terminate();
  mgutils.clear();



  return 0;
  }


MEDCoupling::MEDCouplingFieldDouble *GetFieldComponent( MEDCoupling::MEDCouplingFieldDouble* fieldToExtract, int compID){
    
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> Field= MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
    int nodes = fieldToExtract->getMesh()->getNumberOfNodes();
    MEDCoupling::DataArrayDouble * Arr = fieldToExtract->getArray()->deepCopy();
    MEDCoupling::DataArrayDouble * Arr2 = MEDCoupling::DataArrayDouble::New();
    
    Arr2->alloc(nodes, 1);
    for(int i=0; i<nodes; i++)
        Arr2->setIJ(i,0,Arr->getIJ(i,compID));
    
    Field->setMesh(fieldToExtract->getMesh());
    Field->setArray(Arr2);
    
    Arr->decrRef();
    Arr2->decrRef();
    
    return Field->deepCopy();
}
