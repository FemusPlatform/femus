#include "IbUtils.h"
#include <sstream>
#include <math.h>
#include <iostream>
#include <fstream>

#include "MGFE_conf.h"
#include "MGFE.h"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingMemArray.hxx"



using namespace std;

IbUtils::IbUtils ():InterfaceProjection (), 
_Proc (0),
_ExtractedMesh (NULL),
_ReconInterface (NULL),
_TriangulatedNearBoundary (NULL),
_SolidCellColor (NULL),
_SolidNodeColor (NULL),
_SolidVelocity (NULL),
_FluidCellColor (NULL),
_FluidNodeColor (NULL),
_ExtractedColor (NULL),
_ExtractedVelocity (NULL),
_InterpolatedSolidVelocity (NULL),
_TriangulatedNodesIds (NULL), _CellToExtract (NULL), _InterfaceNormal (NULL)
{
  _SolidBodyMesh = NULL;
  _FluidMesh = NULL;
  read_file ();
}

IbUtils::IbUtils (MEDCoupling::MEDCouplingUMesh * SolidBodyMesh,
		  const MEDCoupling::MEDCouplingUMesh * FluidMesh, int proc):
InterfaceProjection (SolidBodyMesh, FluidMesh, proc, Volume),
_Proc (proc),
_ExtractedMesh (NULL),
_ReconInterface (NULL),
_TriangulatedNearBoundary (NULL),
_SolidCellColor (NULL),
_SolidNodeColor (NULL),
_SolidVelocity (NULL),
_FluidCellColor (NULL),
_FluidNodeColor (NULL),
_ExtractedColor (NULL),
_ExtractedVelocity (NULL),
_InterpolatedSolidVelocity (NULL),
_TriangulatedNodesIds (NULL),
_CellToExtract (NULL),
_InterfaceNormal (NULL)
{
  int dim = FluidMesh->getSpaceDimension ();

  _PieceWiseInterpScheme = INTERP_KERNEL::Triangulation;
  read_file ();
  SetMeshes (SolidBodyMesh, FluidMesh);
  for (int i = 0; i < DIMENSION; i++)
    _Velocity[i] = 0.;
  _iter = 0;
}

IbUtils::IbUtils (MEDCoupling::MEDCouplingUMesh * SolidBodyMesh,
		  const MEDCoupling::MEDCouplingUMesh * FluidMesh, int p_femus, int proc):
InterfaceProjection (SolidBodyMesh, FluidMesh, proc, Volume),
_Proc (proc),
_ExtractedMesh (NULL),
_ReconInterface (NULL),
_TriangulatedNearBoundary (NULL),
_SolidCellColor (NULL),
_SolidNodeColor (NULL),
_SolidVelocity (NULL),
_FluidCellColor (NULL),
_FluidNodeColor (NULL),
_ExtractedColor (NULL),
_ExtractedVelocity (NULL),
_InterpolatedSolidVelocity (NULL),
_TriangulatedNodesIds (NULL),
_CellToExtract (NULL),
_InterfaceNormal (NULL)
{
  int dim = FluidMesh->getSpaceDimension ();

  _PieceWiseInterpScheme = INTERP_KERNEL::Triangulation;
  read_file (p_femus);
  SetMeshes (SolidBodyMesh, FluidMesh);
  for (int i = 0; i < DIMENSION; i++)
    _Velocity[i] = 0.;
  _iter = 0;
}

IbUtils::~IbUtils ()
{

}


void
IbUtils::SetMeshes (MEDCoupling::MEDCouplingUMesh * SolidBodyMesh,
		    const MEDCoupling::MEDCouplingUMesh * FluidMesh)
{
  _SolidBodyMesh = SolidBodyMesh->deepCopy ();
  _FluidMesh = FluidMesh->deepCopy ();

  MEDCoupling::MEDCouplingUMesh * FemusMesh2 = FluidMesh->deepCopy ();
  FemusMesh2->convertQuadraticCellsToLinear ();
  _LinearFluidMesh = FemusMesh2->deepCopy ();

  MEDCoupling::MEDCouplingUMesh * PalaMesh = SolidBodyMesh->deepCopy ();
  PalaMesh->convertQuadraticCellsToLinear ();	// palamesh
  _LinearSolidBodyMesh = PalaMesh->deepCopy ();
  
  
  _LinearFluidMesh->zipCoords();
  _LinearSolidBodyMesh->zipCoords();
  
//   if(_Proc==0){
//   MEDCoupling::WriteUMesh("linear_solid.med",_LinearSolidBodyMesh,true);
//   MEDCoupling::WriteUMesh("linear_fluid.med",_LinearFluidMesh,true);
//   }
//   PrintMed(_LinearFluidMesh    ,"linear_fluid.med");

  FemusMesh2->decrRef ();
  PalaMesh->decrRef ();

  InitColor ();

  MEDCoupling::MEDCouplingRemapper Remap;
  Remap.setPrecision (1e-11);
  Remap.setIntersectionType (_PieceWiseInterpScheme);	// <- more accurate than triangulate
  Remap.prepare (_LinearSolidBodyMesh, _LinearFluidMesh, "P0P0");

//   PrintMed (_SolidNodeColor, "Solid.med", true);

  _FluidNodeColor = InterpolatedField (_SolidNodeColor, 0., 2);
  _FluidCellColor = Remap.transferField (_SolidCellColor, 0.);
  _RemapMatrix.clear ();
  _RemapMatrix = Remap.getCrudeMatrix ();


  if (_Proc == 0)
    MEDCoupling::WriteField ("cell_color.med", _FluidCellColor, true);

  _SolidVelocity =
    MEDCoupling::MEDCouplingFieldDouble::New (MEDCoupling::ON_NODES);
  MEDCoupling::DataArrayDouble * VelArray =
    MEDCoupling::DataArrayDouble::New ();
  VelArray->alloc (_SolidBodyMesh->getNumberOfNodes (), DIMENSION);
  VelArray->fillWithZero ();
  _SolidVelocity->setArray (VelArray);
  _SolidVelocity->setMesh (_SolidBodyMesh);
  _SolidVelocity->setName ("SolidVelocity");
  _SolidVelocity->setNature (MEDCoupling::IntensiveMaximum);
  _SolidVelocity->getArray ()->setInfoOnComponent (0, "x");
  _SolidVelocity->getArray ()->setInfoOnComponent (1, "y");
  VelArray->decrRef ();


  _TriangulatedNodesIds =
    MEDCoupling::MEDCouplingFieldDouble::New (MEDCoupling::ON_NODES);
  return;
}

void
IbUtils::InitColor ()
{

  MEDCoupling::DataArrayDouble * MuFieldArray2 =
    MEDCoupling::DataArrayDouble::New ();
  MuFieldArray2->alloc (_LinearSolidBodyMesh->getNumberOfCells (), 1);
  MuFieldArray2->fillWithValue (1.);
  _SolidCellColor =
    MEDCoupling::MEDCouplingFieldDouble::New (MEDCoupling::ON_CELLS);
  _SolidCellColor->setArray (MuFieldArray2);
  _SolidCellColor->setMesh (_LinearSolidBodyMesh);
  _SolidCellColor->setName ("Cell_Viscosity");
  _SolidCellColor->setNature (MEDCoupling::IntensiveConservation);

  MuFieldArray2->decrRef ();

  MEDCoupling::DataArrayDouble * MuFieldArray =
    MEDCoupling::DataArrayDouble::New ();
  MuFieldArray->alloc (_SolidBodyMesh->getNumberOfNodes (), 1);
  MuFieldArray->fillWithValue (1.);
  _SolidNodeColor =
    MEDCoupling::MEDCouplingFieldDouble::New (MEDCoupling::ON_NODES);
  _SolidNodeColor->setArray (MuFieldArray);
  _SolidNodeColor->setMesh (_SolidBodyMesh);
  _SolidNodeColor->setName ("Node_Viscosity");

  MuFieldArray->decrRef ();

  return;
}

void
IbUtils::CalcInterface (MEDCoupling::MEDCouplingFieldDouble * ProjectedColor)
{
  int Conn[4][4];
  Conn[1][0] = 4;
  Conn[1][1] = 1;
  Conn[1][2] = 5;
  Conn[1][3] = 8;
  Conn[2][0] = 8;
  Conn[2][1] = 5;
  Conn[2][2] = 2;
  Conn[2][3] = 6;
  Conn[3][0] = 7;
  Conn[3][1] = 8;
  Conn[3][2] = 6;
  Conn[3][3] = 3;
  Conn[0][0] = 0;
  Conn[0][1] = 4;
  Conn[0][2] = 8;
  Conn[0][3] = 7;

  _TriangulatedCellsPerSourceCell.clear ();
  _TriangulatedCellsPerSourceCell.clear ();
  _iter++;
  double *Color =
    const_cast < double *>(_FluidCellColor->getArray ()->getPointer ());
  int cells = _FluidCellColor->getMesh ()->getNumberOfCells ();
  std::vector < int >InterfaceCells;
  _SolidArea = 0.;
  for (int i = 0; i < cells; i++)
    {
      if (Color[i] > 1.1)
	{
	  Color[i] *= 0.5;
	  std::cerr << "IBMover: Cell " << i << "  Color Value  " << Color[i]
	    * 2. << std::endl;
	}
      if (Color[i] > _LowerThreshold && Color[i] < _UpperThreshold)
	InterfaceCells.push_back (i);
      _SolidArea += Color[i] * 0.0008;
    }

  std::cout << "Solid Area " << _SolidArea << std::endl;

  int *CellId = new int[InterfaceCells.size ()];

  delete[]_CellToExtract;
  delete[]_InterfaceNormal;

  _CellToExtract = new int[InterfaceCells.size ()];
  _InterfaceNormal = new double[DIMENSION * InterfaceCells.size ()];
  _NumCellToExtract = InterfaceCells.size ();
  for (int i = 0; i < InterfaceCells.size (); i++)
    {
      CellId[i] = InterfaceCells[i];
      _CellToExtract[i] = InterfaceCells[i];
    }
  if (_ExtractedColor != NULL)
    _ExtractedColor->decrRef ();
  _ExtractedColor =
    ProjectedColor->buildSubPart (CellId, CellId + InterfaceCells.size ());
  double *ExNodeColor =
    const_cast < double *>(_ExtractedColor->getArray ()->getPointer ());
  if (_ExtractedMesh != NULL)
    _ExtractedMesh->decrRef ();
  _ExtractedMesh =
    _FluidMesh->buildPartOfMySelf (CellId, CellId + InterfaceCells.size ());
  _ExtractedMesh->zipCoords ();
  _ExtractedMesh->setName ("InterfaceCells");

  int dim = ProjectedColor->getMesh ()->getSpaceDimension ();

  const int WeightDim = dim - 1;
  const int QPoints = _fe[2]->_NoGauss1[WeightDim];
  const int CNodes = _FluidMesh->getNumberOfNodesInCell (0);

  double *Coordinates = new double[dim * CNodes];
  double *FullCoordinates = new double[dim * (2 + CNodes)];
  double *ColVal = new double[CNodes];
  double *dphi = new double[CNodes * dim];
  double InvJac[DIMENSION * DIMENSION];

  std::vector < double >InterfaceCoords;
  std::vector < double >TriangulationCoords;
  std::vector < int >SourceIDs;
  std::vector < int >NodesPerCell;
  NodesPerCell.push_back (0);
  for (int j = 0; j < InterfaceCells.size (); j++)
    {
      std::vector < int >SourceConnectivity;
      std::vector < int >Connectivity;
      std::vector < double >NodeCoordinates;
      _ExtractedMesh->getNodeIdsOfCell (j, Connectivity);
      _FluidMesh->getNodeIdsOfCell (InterfaceCells[j], SourceConnectivity);
      _CoordMatrix.clear ();
      for (int k = 0; k < CNodes; k++)
	{
	  ColVal[k] = ExNodeColor[Connectivity[k]];
	  _ExtractedMesh->getCoordinatesOfNode (Connectivity[k],
						NodeCoordinates);
	  for (int d = 0; d < dim; d++)
	    {
	      Coordinates[k + d * CNodes] = NodeCoordinates[d];
	      FullCoordinates[k + d * (CNodes + 2)] = NodeCoordinates[d];
	      if (k < 4)
		_CoordMatrix.push_back (NodeCoordinates[d]);
	    }

	  NodeCoordinates.clear ();
	}
      Connectivity.clear ();
      double nx = 0., ny = 0., mod = 0.;
      for (int qp = 0; qp < QPoints; qp++)
	{
	  double det = _fe[2]->Jac (qp, Coordinates, InvJac);	//
	  double JxW_g = det * _fe[2]->_weight1[WeightDim][qp];	// weight
	  _fe[2]->get_dphi_gl_g (DIMENSION, qp, InvJac, dphi);	// global coord deriv
	  double dx_g = 0., dy_g = 0.;
	  for (int n = 0; n < CNodes; n++)
	    {
	      dx_g += JxW_g * ColVal[n] * dphi[n];
	      dy_g += JxW_g * ColVal[n] * dphi[n + CNodes];
	    }
	  nx += JxW_g * dx_g;
	  ny += JxW_g * dy_g;
	}
      double norm[DIMENSION];
      mod = sqrt (nx * nx + ny * ny);
      norm[0] = -nx / (mod);
      norm[1] = -ny / (mod);
      _InterfaceNormal[j * DIMENSION] = norm[0];
      _InterfaceNormal[j * DIMENSION + 1] = norm[1];

//      std::cout<<norm[0]*norm[0] + norm[1]*norm[1]<<"  "<<mod<<std::endl;

      double bound_box[4][2];
      double xmin =
	(Coordinates[0] < Coordinates[2]) ? Coordinates[0] : Coordinates[2];
      double xmax =
	(Coordinates[0] > Coordinates[2]) ? Coordinates[0] : Coordinates[2];
      double ymin =
	(Coordinates[0 + CNodes] <
	 Coordinates[2 + CNodes]) ? Coordinates[0 + CNodes] : Coordinates[2 +
									  CNodes];
      double ymax =
	(Coordinates[0 + CNodes] >
	 Coordinates[2 + CNodes]) ? Coordinates[0 + CNodes] : Coordinates[2 +
									  CNodes];

      bound_box[0][0] = bound_box[3][0] = xmin;
      bound_box[1][0] = bound_box[2][0] = xmax;
      bound_box[0][1] = bound_box[1][1] = ymin;
      bound_box[2][1] = bound_box[3][1] = ymax;

      double dxm = xmax - xmin;
      double dym = ymax - ymin;
      double asp_ratio = dym / dxm;
      double dir_ny = (norm[1] > 0) ? 1 : -1;
      double dir_nx = (norm[0] > 0) ? 1 : -1;
      double pen = fabs (norm[0]) / (fabs (norm[1]) + 1.e-10);

      const int Case = (pen > asp_ratio) ? 0 : 1;
      const double tri_area =
	(Case == 0) ? dym * dym * 0.5 / pen : dxm * dxm * 0.5 * pen;
      const double rec_area = dxm * dym;
      const double par_area = rec_area - 2. * tri_area;
      const double col_area = Color[InterfaceCells[j]] * rec_area;
      int case2 = (col_area < tri_area) ? 0 : 1;
      if (case2 == 1)
	case2 = (col_area > (par_area + tri_area)) ? 2 : 1;

      double dx1, dx2, dy1, dy2;
      double dx = 0., dy = 0.;

      int fx = (norm[1] > 0 - 1.e-10) ? 0 : 1;
      int fy = (norm[1] * norm[0] > 0 - 1.e-10) ? 0 : 1;
      int node = 2 * fx + fy;

      if (fabs (norm[0] * norm[1]) > 1.e-3)
	{
	  if (case2 == 0)
	    {
	      dy = sqrt (2. * col_area * pen);
	      dx = dy / pen;
	      dx1 = dx;
	      dy1 = 0.;
	      dx2 = 0.;
	      dy2 = dy;
	    }
	  else if (case2 == 1)
	    {
	      double area = col_area - tri_area;
	      if (Case == 0)
		{
		  dx = area / dym;
		  dy = 0.;
		  dx1 = dx + dym / pen;
		  dy1 = dy;
		  dx2 = dx;
		  dy2 = dym;
		}
	      else
		{
		  dy = area / dxm;
		  dx = 0.;
		  dx1 = dxm;
		  dy1 = dy;
		  dx2 = dx;
		  dy2 = dy + dxm * pen;
		}
	    }
	  else
	    {
	      double area = col_area - (par_area + tri_area);
	      double area2 = tri_area - area;
	      dy = sqrt (area2 * 2. * pen);
	      dx = dy / pen;
	      dx1 = dxm;
	      dy1 = dym - dy;
	      dx2 = dxm - dx;
	      dy2 = dym;
	    }
	  dx1 *= dir_nx;
	  dx2 *= dir_nx;
	  dy1 *= dir_ny;
	  dy2 *= dir_ny;
	}
      else
	{			// rectangular area
	  if (fabs (norm[0]) > 0.5)
	    {
	      dx1 = dx2 = dir_nx * col_area / dym;
	      dy1 = 0;
	      dy2 = dym;
	      node = (norm[0] > 0) ? 0 : 1;
	    }
	  else
	    {
	      dx1 = 0.;
	      dx2 = dxm;
	      dy1 = dy2 = dir_ny * col_area / dxm;
	      node = (norm[1] > 0) ? 0 : 3;
	    }
	}

      double point1[DIMENSION];
      double point2[DIMENSION];
      point1[0] = bound_box[node][0] + dx1;
      point1[1] = bound_box[node][1] + dy1;
      point2[0] = bound_box[node][0] + dx2;
      point2[1] = bound_box[node][1] + dy2;
      _FamilyType = 1;
      _SrcCoordInterpNodes = 4;
      XiEtaChiCalc (point1);
      double CanPos[2];
      CanPos[0] = -sqrt (3. / 5.) /*_XiEtaChi[0]*/ ;
      CanPos[1] = -sqrt (3. / 5.) /*_XiEtaChi[1]*/ ;
      double dphi[18];
      _fe[2]->get_dphi_on_given_node (DIMENSION, Coordinates, CanPos, dphi);
      _XiEtaChi.clear ();
      XiEtaChiCalc (point2);

      InterfaceCoords.push_back (point1[0]);
      InterfaceCoords.push_back (point2[0]);
      InterfaceCoords.push_back (point1[1]);
      InterfaceCoords.push_back (point2[1]);
      FullCoordinates[NDOF_FEM] = point1[0];
      FullCoordinates[NDOF_FEM + 1] = point2[0];
      FullCoordinates[NDOF_FEM + (NDOF_FEM + 2)] = point1[1];
      FullCoordinates[NDOF_FEM + 1 + (NDOF_FEM + 2)] = point2[1];


      // calculation of active nodes
      std::vector < int >fluid_nodes;
      std::vector < int >solid_nodes;
      int center_included = 0;
      int Order[8] = { 0, 4, 1, 5, 2, 6, 3, 7 };
      int Reverse[8] = { 0, 2, 5, 7, 1, 4, 6, 8 };
      int Flag[8];
      std::vector < int >OrderedNodes;
      int LastNo = -1;
      for (int l = 0; l < NDOF_FEM; l++)
	{			// CALCULATION OF FLUID NODES: scalar product between distance node-interaface node and interface normal
	  double scal_prod = 0.;
	  for (int dir = 0; dir < DIMENSION; dir++)
	    {
	      scal_prod +=
		(Coordinates[l + dir * NDOF_FEM] - point1[dir]) * norm[dir];
	    }
	  if (scal_prod > 0.)
	    {
	      fluid_nodes.push_back (l);
	      if (l == 8)
		center_included = 1;
//                 for ( int dir=0; dir<DIMENSION; dir++ ) TriangulationCoords.push_back ( Coordinates[l +  dir*NDOF_FEM] );
	      if (l != 8)
		Flag[l] = 1;
	    }
	  else
	    {
	      solid_nodes.push_back (l);
	      _SolidIds[InterfaceCells[j]].push_back (SourceConnectivity[l]);
	      if (l != 8)
		{
		  Flag[l] = 0;
		  LastNo = (Reverse[l] > LastNo) ? Reverse[l] : LastNo;
		}
	    }
	}

      for (int l = LastNo; l < NDOF_FEM - 1 + LastNo; l++)
	{			// ordering elem connectivity from node close to interface
	  int pos = l % (NDOF_FEM - 1);
	  if (Flag[Order[pos]] == 1)
	    OrderedNodes.push_back (Order[pos]);
	}

      std::vector < int >FinalConn;
      int rectangular_subdomains, act_fluid_nodes;
      if (center_included == 1)
	{			// SEARCHING FOR QUAD4 SUB ELEMENTS
	  std::vector < int >fluid_quad4;
	  for (int rec = 0; rec < 4; rec++)
	    {
	      int fluid_rec = 1;
	      for (int node = 0; node < 4; node++)
		{
		  for (int r = 0; r < solid_nodes.size (); r++)
		    {
		      if (Conn[rec][node] == solid_nodes[r])
			{
			  fluid_rec = 0;
			  break;
			}
		    }
		}
	      if (fluid_rec == 1)
		fluid_quad4.push_back (rec);
	    }
	  if (fluid_quad4.size () > 0)
	    {
	      int numRec = fluid_quad4.size ();
	      int FirstPos = -1, LastPos = -1;
	      int Nodes = OrderedNodes.size ();
	      for (int n = 0; n < Nodes; n++)
		{
		  int RecNode = 1;
		  for (int r = 0; r < numRec; r++)
		    {
		      if (fluid_quad4[r] == OrderedNodes[n])
			{
			  FirstPos = (FirstPos < 0) ? n : FirstPos;
			  LastPos = (n >= FirstPos) ? n : LastPos;
			}
		    }
		}
	      std::vector < int >ElemenConn = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
	      for (int r = 0; r < numRec; r++)
		{		// ADDING RECTANGULAR SUB ELEMENTS
		  std::vector < int >Sequence;
		  for (int node = 0; node < 4; node++)
		    {
		      Sequence.push_back (Conn[fluid_quad4[r]][node]);
		      SourceIDs.push_back (SourceConnectivity
					   [Conn[fluid_quad4[r]][node]]);
		    }
		  InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
					  ElemenConn, Sequence,
					  TriangulationCoords, Coordinates,
					  NDOF_FEM);
//                     _TriangulatedCellsPerSourceCell[InterfaceCells[j]].push_back ( NodesPerCell.size() );
//                     NodesPerCell.push_back ( NodesPerCell[NodesPerCell.size() -1] + 4 );
//                     for ( int node=0; node<4; node++ )
//                         for ( int dim=0; dim<DIMENSION; dim++ )
//                             TriangulationCoords.push_back ( Coordinates[Conn[fluid_quad4[r]][node] +  dim*NDOF_FEM] );
		}
	      // CUTTING CONNECTIVITY VECTOR -> ELIMINATING QUAD4 CORNERS AND ADDITION OF ELEMENT MID POINT
	      for (int n = 0; n < Nodes; n++)
		{
		  if (n < FirstPos)
		    FinalConn.push_back (OrderedNodes[n]);
		  else if (n == FirstPos)
		    FinalConn.push_back (8);
		  else if (n > LastPos)
		    FinalConn.push_back (OrderedNodes[n]);
		}
	    }
	  rectangular_subdomains = fluid_quad4.size ();
	  fluid_quad4.clear ();
	}
      else
	{			// no presence of rectangles
	  rectangular_subdomains = 0;
	  for (int l = 0; l < OrderedNodes.size (); l++)
	    FinalConn.push_back (OrderedNodes[l]);
	}

      OrderedNodes.clear ();
      OrderedNodes.resize (FinalConn.size () + 2);
      int First, Last;
      if (fabs (FullCoordinates[FinalConn[0]] - FullCoordinates[NDOF_FEM]) <
	  1.e-10
	  || fabs (FullCoordinates[FinalConn[0] + NDOF_FEM + 2] -
		   FullCoordinates[NDOF_FEM + NDOF_FEM + 2]) < 1.e-10)
	{
	  First = NDOF_FEM;
	  Last = NDOF_FEM + 1;
	  OrderedNodes[0] = First;
	  OrderedNodes[FinalConn.size () + 1] = Last;
	}
      else
	{
	  First = NDOF_FEM + 1;
	  Last = NDOF_FEM;
	  OrderedNodes[0] = First;
	  OrderedNodes[FinalConn.size () + 1] = Last;
	}
      for (int i = 0; i < FinalConn.size (); i++)
	OrderedNodes[i + 1] = FinalConn[i];


      act_fluid_nodes = FinalConn.size ();
//         std::cout<<"Polygonal cell with "<<rectangular_subdomains<<" linear quad4 and "<<act_fluid_nodes<<" remaining fluid nodes \n";
//         for ( int i=0; i<FinalConn.size() +2; i++ ) std::cout<<OrderedNodes[i]<<" ";
//         std::cout<<"\n";
//         // ADDING INTERFACE NODE CLOSEST TO CONNECTIVITY FIRST POINT
//         if ( fabs ( Coordinates[FinalConn[0]] - point1[0] ) <1.e-10 ||
//                 fabs ( Coordinates[FinalConn[0]+NDOF_FEM] - point1[1] ) <1.e-10
//            ) {
//             TriangulationCoords.push_back ( point1[0] );
//             TriangulationCoords.push_back ( point1[1] );
//         } else {
//             TriangulationCoords.push_back ( point2[0] );
//             TriangulationCoords.push_back ( point2[1] );
//         }
//         // ADDING CONNECTIVITY NODES
//         for ( int k=0; k<FinalConn.size(); k++ )
//             for ( int dir=0; dir<DIMENSION; dir++ ) {
//                 TriangulationCoords.push_back ( Coordinates[FinalConn[k] +  dir*NDOF_FEM] );
//             }
//         // ADDING INTERFACE NODE CLOSEST TO CONNECTIVITY LAST POINT
//         if ( fabs ( Coordinates[FinalConn[FinalConn.size()-1]] - point1[0] ) <1.e-10 ||
//                 fabs ( Coordinates[FinalConn[FinalConn.size()-1]+NDOF_FEM] - point1[1] ) <1.e-10
//            ) {
//             TriangulationCoords.push_back ( point1[0] );
//             TriangulationCoords.push_back ( point1[1] );
//         } else {
//             TriangulationCoords.push_back ( point2[0] );
//             TriangulationCoords.push_back ( point2[1] );
//         }
//         // NUMBER OF NODES
//         NodesPerCell.push_back ( NodesPerCell[NodesPerCell.size() -1] + FinalConn.size() + 2 );

// // //         SPLITTING REMAINING NODES
      int offset = NDOF_FEM + 2;
      switch (act_fluid_nodes)
	{
	case 5:
	  {
	    std::vector < int >Sequence1
	    {
	    0, 1, 2, 3};
	    std::vector < int >Sequence2
	    {
	    3, 4, 5, 6};
	    std::vector < int >Sequence3
	    {
	    0, 3, 6};
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence1,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence2,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence3,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    for (int j = 0; j < Sequence1.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence1[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence1[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence1[j]]]);
	    for (int j = 0; j < Sequence2.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence2[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence2[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence2[j]]]);
	    for (int j = 0; j < Sequence3.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence3[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence3[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence3[j]]]);
	    Sequence1.clear ();
	    Sequence2.clear ();
	    Sequence3.clear ();
	  }
	  break;
	case 4:
	  {
	    std::vector < int >Sequence1, Sequence2, Sequence3;
	    if (rectangular_subdomains == 2)
	      {
		Sequence1.resize (3);
		Sequence1 =
		{
		0, 1, 2};
		Sequence2.resize (4);
		Sequence2 =
		{
		0, 2, 3, 5};
		Sequence3.resize (3);
		Sequence3 =
		{
		3, 4, 5};
	      }
	    else if (rectangular_subdomains == 1)
	      {
		if (OrderedNodes[2] == 8)
		  {
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    0, 1, 2};
		    Sequence2.resize (4);
		    Sequence2 =
		    {
		    2, 3, 4, 5};
		    Sequence3.resize (3);
		    Sequence3 =
		    {
		    0, 2, 5};
		  }
		else if (OrderedNodes[3] == 8)
		  {
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    3, 4, 5};
		    Sequence2.resize (4);
		    Sequence2 =
		    {
		    0, 1, 2, 3};
		    Sequence3.resize (3);
		    Sequence3 =
		    {
		    0, 3, 5};
		  }
		else
		  {
		    printf
		      ("IButils::Interface error: unknown combination with %d free nodes and %d rectangles\n",
		       act_fluid_nodes, rectangular_subdomains);
		  }
	      }
	    else if (rectangular_subdomains == 0)
	      {
		if (fabs
		    (FullCoordinates[OrderedNodes[0]] -
		     FullCoordinates[OrderedNodes[2]]) < 1.e-10
		    || fabs (FullCoordinates[OrderedNodes[0] + offset] -
			     FullCoordinates[OrderedNodes[2] + offset]) <
		    1.e-10)
		  {
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    3, 4, 5};
		    Sequence2.resize (4);
		    Sequence2 =
		    {
		    0, 1, 3, 5};
		    Sequence3.resize (3);
		    Sequence3 =
		    {
		    1, 2, 3};
		  }
		else
		  {
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    0, 1, 2};
		    Sequence2.resize (4);
		    Sequence2 =
		    {
		    0, 2, 4, 5};
		    Sequence3.resize (3);
		    Sequence3 =
		    {
		    2, 3, 4};
		  }
	      }
	    else
	      {
		printf
		  ("IButils::Interface error: unknown combination with %d free nodes and %d rectangles\n",
		   act_fluid_nodes, rectangular_subdomains);
	      }
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence1,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence2,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence3,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    for (int j = 0; j < Sequence1.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence1[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence1[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence1[j]]]);
	    for (int j = 0; j < Sequence2.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence2[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence2[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence2[j]]]);
	    for (int j = 0; j < Sequence3.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence3[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence3[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence3[j]]]);

	    Sequence1.clear ();
	    Sequence2.clear ();
	    Sequence3.clear ();
	  }
	  break;
	case 3:
	  {
	    std::vector < int >Sequence1, Sequence2;
	    if (rectangular_subdomains == 0)
	      {
		if (fabs
		    (FullCoordinates[OrderedNodes[1]] -
		     FullCoordinates[OrderedNodes[3]]) < 1.e-10
		    || fabs (FullCoordinates[OrderedNodes[1] + offset] -
			     FullCoordinates[OrderedNodes[3] + offset]) <
		    1.e-10)
		  {
		    std::vector < int >Sequence3 = { 0, 1, 2 };
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    0, 2, 4};
		    Sequence2.resize (3);
		    Sequence2 =
		    {
		    4, 2, 3};
		    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
					    OrderedNodes, Sequence3,
					    TriangulationCoords,
					    FullCoordinates, offset);
		    for (int j = 0; j < Sequence3.size (); j++)
		      SourceIDs.push_back ((OrderedNodes[Sequence3[j]] >=
					    NDOF_FEM) ? -(1 +
							  OrderedNodes
							  [Sequence3[j]] %
							  NDOF_FEM) :
					   SourceConnectivity[OrderedNodes
							      [Sequence3
							       [j]]]);
		    Sequence3.clear ();
		  }
		else
		  {
		    Sequence1.resize (3);
		    Sequence1 =
		    {
		    1, 3, 2};
		    Sequence2.resize (4);
		    Sequence2 =
		    {
		    0, 4, 3, 1};
		  }
	      }
	    else if (rectangular_subdomains == 2
		     || rectangular_subdomains == 3)
	      {
		std::vector < int >Sequence3 = { 0, 1, 2 };
		Sequence1.resize (3);
		Sequence1 =
		{
		0, 2, 4};
		Sequence2.resize (3);
		Sequence2 =
		{
		4, 2, 3};
		InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
					OrderedNodes, Sequence3,
					TriangulationCoords, FullCoordinates,
					offset);
		for (int j = 0; j < Sequence3.size (); j++)
		  SourceIDs.push_back ((OrderedNodes[Sequence3[j]] >=
					NDOF_FEM) ? -(1 +
						      OrderedNodes[Sequence3
								   [j]] %
						      NDOF_FEM) :
				       SourceConnectivity[OrderedNodes
							  [Sequence3[j]]]);
		Sequence3.clear ();
	      }
	    else
	      {
		printf
		  ("IButils::Interface error: unknown combination with %d free nodes and %d rectangles\n",
		   act_fluid_nodes, rectangular_subdomains);
	      }
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence1,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence2,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    for (int j = 0; j < Sequence1.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence1[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence1[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence1[j]]]);
	    for (int j = 0; j < Sequence2.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence2[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence2[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence2[j]]]);
	    Sequence1.clear ();
	    Sequence2.clear ();
	  }
	  break;
	case 2:
	  {
	    // check on connectivity
	    if ((fabs
		 (FullCoordinates[OrderedNodes[0]] -
		  FullCoordinates[OrderedNodes[1]]) > 1.e-10
		 && fabs (FullCoordinates[OrderedNodes[0] + offset] -
			  FullCoordinates[OrderedNodes[1] + offset]) > 1.e-10)
		||
		(fabs
		 (FullCoordinates[OrderedNodes[2]] -
		  FullCoordinates[OrderedNodes[3]]) > 1.e-10
		 && fabs (FullCoordinates[OrderedNodes[2] + offset] -
			  FullCoordinates[OrderedNodes[3] + offset]) >
		 1.e-10))
	      {
		int first = OrderedNodes[0];
		OrderedNodes[0] = OrderedNodes[3];
		OrderedNodes[3] = first;
	      }
	    // find aligned nodes
	    std::vector < int >Sequence1 (3), Sequence2 (3);
	    if (fabs
		(FullCoordinates[OrderedNodes[0]] -
		 FullCoordinates[OrderedNodes[2]]) < 1.e-10
		|| fabs (FullCoordinates[OrderedNodes[0] + offset] -
			 FullCoordinates[OrderedNodes[2] + offset]) < 1.e-10)
	      {
		Sequence1 =
		{
		0, 1, 3};
		Sequence2 =
		{
		1, 2, 3};
	      }
	    else
	      {
		Sequence1 =
		{
		0, 1, 2};
		Sequence2 =
		{
		0, 2, 3};
	      }
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence1,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence2,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    for (int j = 0; j < Sequence1.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence1[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence1[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence1[j]]]);
	    for (int j = 0; j < Sequence2.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence2[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence2[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence2[j]]]);
	    Sequence1.clear ();
	    Sequence2.clear ();

	  }
	  break;
	case 1:
	  {
	    std::vector < int >Sequence1 = { 0, 1, 2 };
	    InsertTriangulatedCell (InterfaceCells[j], NodesPerCell,
				    OrderedNodes, Sequence1,
				    TriangulationCoords, FullCoordinates,
				    offset);
	    for (int j = 0; j < Sequence1.size (); j++)
	      SourceIDs.push_back ((OrderedNodes[Sequence1[j]] >=
				    NDOF_FEM) ? -(1 +
						  OrderedNodes[Sequence1[j]] %
						  NDOF_FEM) :
				   SourceConnectivity[OrderedNodes
						      [Sequence1[j]]]);
	    Sequence1.clear ();
	  }
	  break;
	default:
	  std::cout << "Unknown case with " << act_fluid_nodes <<
	    " active fluid nodes \n";
//        int a=1;
	  break;
	}
    }

  int NewCells = NodesPerCell.size () - 1;
  int NewNodes = NodesPerCell[NewCells];
  std::cout << "IBUTILS ncells " << NewCells << std::endl;
  if (_TriangulatedNearBoundary != NULL)
    _TriangulatedNearBoundary->decrRef ();
  _TriangulatedNearBoundary =
    MEDCoupling::MEDCouplingUMesh::New ("RefinedMesh2", 2);
  MEDCoupling::DataArrayDouble * TrianCoordArray =
    MEDCoupling::DataArrayDouble::New ();
  MEDCoupling::DataArrayDouble * SourceConnArray =
    MEDCoupling::DataArrayDouble::New ();
  TrianCoordArray->alloc (NewNodes, DIMENSION);
  SourceConnArray->alloc (NewNodes, 1);
  _TriangulatedNearBoundary->allocateCells (NewCells);
  for (int i = 0; i < NewCells; i++)
    {				// REFINED CELLS ================================
      int nodes = NodesPerCell[i + 1] - NodesPerCell[i];
      int *connect = new int[nodes];
      for (int k = 0; k < nodes; k++)
	{
	  int new_node_id = k + NodesPerCell[i];
	  connect[k] = new_node_id;
	  for (int d = 0; d < DIMENSION; d++)
	    TrianCoordArray->setIJ (new_node_id, d,
				    TriangulationCoords[new_node_id *
							DIMENSION + d]);
	  SourceConnArray->setIJ (new_node_id, 0, SourceIDs[new_node_id]);
	}
      if (nodes == 4)
	_TriangulatedNearBoundary->insertNextCell (INTERP_KERNEL::NORM_QUAD4,
						   nodes, connect);
      else if (nodes == 3)
	_TriangulatedNearBoundary->insertNextCell (INTERP_KERNEL::NORM_TRI3,
						   nodes, connect);
      else
	_TriangulatedNearBoundary->insertNextCell (INTERP_KERNEL::
						   NORM_POLYGON, nodes,
						   connect);
      delete[]connect;
    }
  _TriangulatedNearBoundary->finishInsertingCells ();
  _TriangulatedNearBoundary->setCoords (TrianCoordArray);
//     MEDCoupling::DataArrayInt  *Simplexity = Triangulation->simplexize ( 0 );
//     MEDLoader::WriteUMesh ( "RESU_MED/Triangulation_"+to_string ( _iter ) +".med", _TriangulatedNearBoundary, true );

  _TriangulatedNodesIds->setArray (SourceConnArray);
  _TriangulatedNodesIds->setMesh (_TriangulatedNearBoundary);
  _TriangulatedNodesIds->setName ("NodeIds");

  _TriangulatedNearBoundary->zipCoords ();
//      _TriangulatedNodesIds->getMesh()->zipCoords();
//     if(_Proc==0)MEDLoader::WriteField ( "RESU_MED/Triangulation_"+to_string ( _iter ) +".med", _TriangulatedNodesIds, true );

  int connect[2];
  int Intcells = InterfaceCoords.size () / (DIMENSION * 2);
  int Intnodes = InterfaceCoords.size () / (DIMENSION);

  if (_ReconInterface != NULL)
    _ReconInterface->decrRef ();
  _ReconInterface = MEDCoupling::MEDCouplingUMesh::New ("RefinedMesh", 1);
  MEDCoupling::DataArrayDouble * CoordArray =
    MEDCoupling::DataArrayDouble::New ();
  CoordArray->alloc (Intnodes, DIMENSION);
  _ReconInterface->allocateCells (Intcells);

  for (int i = 0; i < Intcells; i++)
    {				// REFINED CELLS ================================
      for (int k = 0; k < 2; k++)
	{
	  int new_node_id = k + i * 2;
	  connect[k] = new_node_id;
	  for (int d = 0; d < DIMENSION; d++)
	    {
	      CoordArray->setIJ (connect[k], d,
				 InterfaceCoords[2 * i * DIMENSION + k +
						 d * 2]);
	    }
	}
      _ReconInterface->insertNextCell (INTERP_KERNEL::NORM_SEG2, 2, connect);
    }				// END ADDING REFINED CELLS ==============================================================
  _ReconInterface->finishInsertingCells ();
  _ReconInterface->setCoords (CoordArray);
  _ReconInterface->zipCoords ();

  TrianCoordArray->decrRef ();
  SourceConnArray->decrRef ();
  CoordArray->decrRef ();

  delete[]CellId;
  delete[]dphi;
  delete[]Coordinates;
  delete[]ColVal;
  delete[]FullCoordinates;
  return;
}

void
IbUtils::CalcInterface2 (MEDCoupling::MEDCouplingFieldDouble * ProjectedColor)
{
  int Conn[4][4];
  Conn[1][0] = 4;
  Conn[1][1] = 1;
  Conn[1][2] = 5;
  Conn[1][3] = 8;
  Conn[2][0] = 8;
  Conn[2][1] = 5;
  Conn[2][2] = 2;
  Conn[2][3] = 6;
  Conn[3][0] = 7;
  Conn[3][1] = 8;
  Conn[3][2] = 6;
  Conn[3][3] = 3;
  Conn[0][0] = 0;
  Conn[0][1] = 4;
  Conn[0][2] = 8;
  Conn[0][3] = 7;

  _TriangulatedCellsPerSourceCell.clear ();
  _TriangulatedCellsPerSourceCell.clear ();
  _iter++;
  double *Color =
    const_cast < double *>(_FluidCellColor->getArray ()->getPointer ());
  int cells = _FluidCellColor->getMesh ()->getNumberOfCells ();
  std::vector < int >InterfaceCells;
  _SolidArea = 0.;
  for (int i = 0; i < cells; i++)
    {
      if (Color[i] > 1.1)
	{
	  Color[i] *= 0.5;
	  std::cerr << "IBMover: Cell " << i << "  Color Value  " << Color[i]
	    * 2. << std::endl;
	}
      if (Color[i] > _LowerThreshold && Color[i] < _UpperThreshold)
	InterfaceCells.push_back (i);
      _SolidArea += Color[i] * 0.0008;
    }

  std::cout << "Solid Area " << _SolidArea << std::endl;

  int *CellId = new int[InterfaceCells.size ()];

  delete[]_CellToExtract;
  delete[]_InterfaceNormal;

  _CellToExtract = new int[InterfaceCells.size ()];
  _InterfaceNormal = new double[DIMENSION * InterfaceCells.size ()];
  _NumCellToExtract = InterfaceCells.size ();
  for (int i = 0; i < InterfaceCells.size (); i++)
    {
      CellId[i] = InterfaceCells[i];
      _CellToExtract[i] = InterfaceCells[i];
    }
  if (_ExtractedColor != NULL)
    _ExtractedColor->decrRef ();
  _ExtractedColor =
    ProjectedColor->buildSubPart (CellId, CellId + InterfaceCells.size ());
  double *ExNodeColor =
    const_cast < double *>(_ExtractedColor->getArray ()->getPointer ());
  if (_ExtractedMesh != NULL)
    _ExtractedMesh->decrRef ();
  _ExtractedMesh =
    _FluidMesh->buildPartOfMySelf (CellId, CellId + InterfaceCells.size ());
  _ExtractedMesh->zipCoords ();
  _ExtractedMesh->setName ("InterfaceCells");

  int dim = ProjectedColor->getMesh ()->getSpaceDimension ();

  const int WeightDim = dim - 1;
  const int QPoints = _fe[2]->_NoGauss1[WeightDim];
  const int CNodes = _FluidMesh->getNumberOfNodesInCell (0);

  double *Coordinates = new double[dim * CNodes];
  double *FullCoordinates = new double[dim * (2 + CNodes)];
  double *FullCoordinates_relative = new double[dim * (2 + CNodes)];
  double *ColVal = new double[CNodes];
  double *dphi = new double[CNodes * dim];
  double InvJac[DIMENSION * DIMENSION];
  double LinCoord[8];

  std::vector < double >InterfaceCoords;
  std::vector < double >TriangulationCoords;
  std::vector < int >SourceIDs;
  std::vector < int >NodesPerCell;
  NodesPerCell.push_back (0);
  for (int j = 0; j < InterfaceCells.size (); j++)
    {
      std::vector < int >SourceConnectivity;
      std::vector < int >Connectivity;
      std::vector < double >NodeCoordinates;
      _ExtractedMesh->getNodeIdsOfCell (j, Connectivity);
      _FluidMesh->getNodeIdsOfCell (InterfaceCells[j], SourceConnectivity);
      _CoordMatrix.clear ();
      for (int k = 0; k < CNodes; k++)
	{
	  ColVal[k] = ExNodeColor[Connectivity[k]];
	  _ExtractedMesh->getCoordinatesOfNode (Connectivity[k],
						NodeCoordinates);
	  for (int d = 0; d < dim; d++)
	    {
	      Coordinates[k + d * CNodes] = NodeCoordinates[d];
	      FullCoordinates[k + d * (CNodes + 2)] = NodeCoordinates[d];
	      if (k < 4)
		{
		  _CoordMatrix.push_back (NodeCoordinates[d]);
		  LinCoord[k] = NodeCoordinates[0];
		  LinCoord[k + 4] = NodeCoordinates[1];
		}
	    }

	  NodeCoordinates.clear ();
	}
      Connectivity.clear ();
      double nx = 0., ny = 0., mod = 0.;

      // now node 0 is in coordinates (0,0)
      for (int k = 0; k < CNodes; k++)
	{
	  for (int d = 0; d < dim; d++)
	    {
	      FullCoordinates_relative[k + d * (CNodes + 2)] =
		FullCoordinates[k + d * (CNodes + 2)] -
		FullCoordinates[d * (CNodes + 2)];
	    }
	}

      double Area = 0;
      // CALCULATION OF NORMAL VECTOR - MEAN VOLUME-FRACTION GRADIENT
      for (int qp = 0; qp < QPoints; qp++)
	{
	  double det = _fe[2]->Jac (qp, Coordinates, InvJac);	//
	  double JxW_g = det * _fe[2]->_weight1[WeightDim][qp];	// weight
	  _fe[2]->get_dphi_gl_g (DIMENSION, qp, InvJac, dphi);	// global coord deriv
	  double dx_g = 0., dy_g = 0.;
	  for (int n = 0; n < CNodes; n++)
	    {
	      dx_g += JxW_g * ColVal[n] * dphi[n];
	      dy_g += JxW_g * ColVal[n] * dphi[n + CNodes];
	    }
	  nx += JxW_g * dx_g;
	  ny += JxW_g * dy_g;
	  Area += JxW_g;
	}
      double norm[DIMENSION], new_norm[DIMENSION], int_incl;
      mod = sqrt (nx * nx + ny * ny);
      norm[0] = -nx / (mod);
      norm[1] = -ny / (mod);

      double Area2 = QuadrangleArea (LinCoord);
      Area = Area2;

      _InterfaceNormal[j * DIMENSION] = norm[0];
      _InterfaceNormal[j * DIMENSION + 1] = norm[1];

      int_incl = -norm[0] / norm[1];

      // angle between 0 and 1
      double dx01 = FullCoordinates_relative[1];
      double dy01 = FullCoordinates_relative[1 + (CNodes + 2)];
      double r01 = sqrt (dy01 * dy01 + dx01 * dx01);
      double cos01 = dx01 / (r01);
      double sin01 = dy01 / (r01);

      // angle between 0 and 1
      double dx03 = FullCoordinates_relative[3];
      double dy03 = FullCoordinates_relative[3 + (CNodes + 2)];
      double r03 = sqrt (dy03 * dy03 + dx03 * dx03);
      double cos03 = dx03 / (r03);
      double sin03 = dy03 / (r03);

//      x' = x cos f - y sin f
//      y' = y cos f + x sin f

      double diag02, diag13, min_diag, max_diag, new_pen;
      diag02 =
	FullCoordinates_relative[2 +
				 (CNodes +
				  2)] / (FullCoordinates_relative[2]);
      diag13 =
	(FullCoordinates_relative[3 + (CNodes + 2)] -
	 FullCoordinates_relative[1 +
				  (CNodes +
				   2)]) / (FullCoordinates_relative[3] -
					   FullCoordinates_relative[1]);

      if (diag02 < diag13)
	{
	  min_diag = diag02;
	  max_diag = diag13;
	}
      else
	{
	  min_diag = diag13;
	  max_diag = diag02;
	}

      int RotCase = 0;
      if (int_incl > min_diag && int_incl < max_diag)
	RotCase = 1;

//         std::cout<<"==========================================================\n";
      // rotation of element

      double min_x, max_x;
      double x2, x3, y2, y3;
      if (RotCase == 0)
	{			// point 1 will be at (1,0) with normalized coordinates
	  for (int k = 0; k < CNodes; k++)
	    {
	      double oldx = FullCoordinates_relative[k];
	      double oldy = FullCoordinates_relative[k + (CNodes + 2)];
	      FullCoordinates_relative[k] = (oldx * cos01 + oldy * sin01);
	      FullCoordinates_relative[k + (CNodes + 2)] =
		(oldy * cos01 - oldx * sin01);
	    }
	  new_norm[0] = norm[0] * cos01 + norm[1] * sin01;
	  new_norm[1] = norm[1] * cos01 - norm[0] * sin01;
	  x2 = FullCoordinates_relative[2];
	  x3 = FullCoordinates_relative[3];
	  y2 = FullCoordinates_relative[2 + (CNodes + 2)];
	  y3 = FullCoordinates_relative[3 + (CNodes + 2)];
	  min_x = 0.;
	  max_x = FullCoordinates_relative[1];
	}
      else
	{			// point 3 will be at (-1,0) with normalized coordinates
	  for (int k = 0; k < CNodes; k++)
	    {
	      double oldx = FullCoordinates_relative[k];
	      double oldy = FullCoordinates_relative[k + (CNodes + 2)];
	      FullCoordinates_relative[k] = (-oldx * cos03 - oldy * sin03);
	      FullCoordinates_relative[k + (CNodes + 2)] =
		(-oldy * cos03 + oldx * sin03);
	    }
	  new_norm[0] = -norm[0] * cos03 - norm[1] * sin03;
	  new_norm[1] = -norm[1] * cos03 + norm[0] * sin03;
	  x2 = FullCoordinates_relative[1];
	  x3 = FullCoordinates_relative[2];
	  y2 = FullCoordinates_relative[1 + (CNodes + 2)];
	  y3 = FullCoordinates_relative[2 + (CNodes + 2)];
	  min_x = FullCoordinates_relative[3];
	  max_x = 0.;
	}
      new_pen = -new_norm[0] / new_norm[1];

      double pen_23 = (y3 - y2) / (x3 - x2);

      double VerT1[6], VerT2[6], Trap[8];
      // vertices for T1
      double tilde_xt1 = -y3 / new_pen + x3;
      double tilde_yt1 = 0;
      Trap[0] = tilde_xt1;
      Trap[4] = 0.;
      Trap[3] = x3;
      Trap[7] = y3;

      int CutTypeT1, CutTypeT2;
      CutTypeT1 = 0;
      if (tilde_xt1 < min_x)
	{
	  tilde_xt1 =
	    (y3 - pen_23 * x3 + new_pen * min_x) / (new_pen - pen_23);
	  tilde_yt1 = new_pen * (tilde_xt1 - min_x);
	  Trap[0] = min_x;
	  Trap[4] = 0.;
	  Trap[3] = tilde_xt1;
	  Trap[7] = tilde_yt1;
	  CutTypeT1 = 1;
	}
      VerT1[0] = x3;
      VerT1[1] = min_x;
      VerT1[2] = tilde_xt1;
      VerT1[3] = y3;
      VerT1[4] = 0;
      VerT1[5] = tilde_yt1;

      // vertices for T2
      double tilde_xt2 = -y2 / new_pen + x2;
      double tilde_yt2 = 0;
      Trap[1] = tilde_xt2;
      Trap[5] = 0.;
      Trap[2] = x2;
      Trap[6] = y2;
      CutTypeT2 = 0;
      if (tilde_xt2 > max_x)
	{
	  tilde_xt2 =
	    (y3 - pen_23 * x3 + new_pen * max_x) / (new_pen - pen_23);
	  tilde_yt2 = new_pen * (tilde_xt2 - max_x);
	  Trap[1] = max_x;
	  Trap[5] = 0.;
	  Trap[2] = tilde_xt2;
	  Trap[6] = tilde_yt2;
	  CutTypeT2 = 1;
	}
      VerT2[0] = x2;
      VerT2[1] = max_x;
      VerT2[2] = tilde_xt2;
      VerT2[3] = y2;
      VerT2[4] = 0;
      VerT2[5] = tilde_yt2;

      double AreaT1 = TriangleArea (VerT1);
      double AreaT2 = TriangleArea (VerT2);
      double TrapArea = Area - AreaT1 - AreaT2;
      const double col_area = Color[InterfaceCells[j]] * Area;

      double dir_ny = (new_norm[1] > 0) ? 1 : -1;
      double dir_nx = (new_norm[0] > 0) ? 1 : -1;

      int AreaCase;
      double net_area;
      if (dir_nx == 1)
	{
	  if (col_area < AreaT1)
	    {			// T1
	      AreaCase = 0;
	      net_area = col_area;
	    }
	  else if (col_area > AreaT1 && col_area < (AreaT1 + TrapArea))
	    {			// Trapezoid
	      AreaCase = 1;
	      net_area = col_area - AreaT1;
	    }
	  else
	    {			// T2
	      AreaCase = 2;
	      net_area = col_area - AreaT1 - TrapArea;
	    }
	}
      else
	{
	  if (col_area < AreaT2)
	    {			// T2
	      AreaCase = 2;
	      net_area = col_area;
	    }
	  else if (col_area > AreaT2 && col_area < (AreaT2 + TrapArea))
	    {			// Trapezoid
	      AreaCase = 1;
	      net_area = col_area - AreaT2;
	    }
	  else
	    {			// T1
	      AreaCase = 0;
	      net_area = col_area - AreaT2 - TrapArea;
	    }
	}

      double intP1[2];
      double intP2[2];

      if (AreaCase == 1)
	{			// Trapezoid
	  double task_area = (dir_nx == 1) ? net_area : TrapArea - net_area;
	  if (fabs (pen_23) > 1.e-8)
	    {
	      double aux_x = -y3 / pen_23 + x3;

	      double SupTriVer[6] =
		{ Trap[0], Trap[3], aux_x, 0, Trap[7], 0 };
	      double AuxArea = TriangleArea (SupTriVer);
	      double sign = (pen_23 > 0) ? 1. : -1.;
	      double dx = Trap[1] - Trap[0];
	      double dy = Trap[6] - Trap[7];

	      double A = dx * dy;
	      double B = Trap[7] * dx + (Trap[0] - aux_x) * dy;
	      double C =
		(Trap[0] - aux_x) * Trap[7] - 2 * (AuxArea +
						   sign * task_area) / sign;

	      double f1 = (-B + sqrt (B * B - 4. * A * C)) / (2. * A);
	      double f2 = (-B - sqrt (B * B - 4. * A * C)) / (2. * A);

	      double f = (fabs (f1) < 1 + 1.e-10) ? f1 : f2;

	      intP1[0] = Trap[0] + dx * f;
	      intP1[1] = 0.;

	      intP2[1] = dy * f + Trap[7];
	      intP2[0] = intP2[1] / new_pen + intP1[0];
	    }
	  else
	    {			// Considered parallelogramm
	      double height = y2;
	      double base = task_area / height;
	      intP1[0] = Trap[0] + base;
	      intP1[1] = 0.;

	      intP2[1] = y2;
	      intP2[0] = intP2[1] / new_pen + intP1[0];
	    }
	}
      if (AreaCase == 0)
	{			// Triangle T1
	  double task_area = (dir_nx == 1) ? net_area : AreaT1 - net_area;

	  double ratio = task_area / AreaT1;

	  if (CutTypeT1 == 0)
	    {			// tildex1 on x0-x1 -- interface crossing P0P3 and P0TildeX
	      double dxP1 = tilde_xt1 - min_x;

	      double dxP2 = x3 - min_x;
	      double dyP2 = y3;

	      intP1[0] = min_x + sqrt (ratio) * dxP1;
	      intP1[1] = 0.;

	      intP2[0] = min_x + sqrt (ratio) * dxP2;
	      intP2[1] = 0. + sqrt (ratio) * dyP2;
	    }
	  if (CutTypeT1 == 1)
	    {
	      double dxP1 = min_x - x3;
	      double dyP1 = 0 - y3;

	      double dxP2 = tilde_xt1 - x3;
	      double dyP2 = tilde_yt1 - y3;

	      intP1[0] = x3 + sqrt (ratio) * dxP1;
	      intP1[1] = y3 + sqrt (ratio) * dyP1;

	      intP2[0] = x3 + sqrt (ratio) * dxP2;
	      intP2[1] = y3 + sqrt (ratio) * dyP2;
	    }
	}
      if (AreaCase == 2)
	{			// Triangle T2
	  double task_area = (dir_nx == 1) ? AreaT2 - net_area : net_area;
	  double ratio = task_area / AreaT2;

	  if (CutTypeT2 == 0)
	    {			// tildex1 on x0-x1 -- interface crossing P1P2 and P1TildeX
	      double dxP1 = tilde_xt2 - max_x;

	      double dxP2 = x2 - max_x;
	      double dyP2 = y2;

	      intP1[0] = max_x + sqrt (ratio) * dxP1;
	      intP1[1] = 0.;

	      intP2[0] = max_x + sqrt (ratio) * dxP2;
	      intP2[1] = 0. + sqrt (ratio) * dyP2;
//                 std::cout<<" cut0, start pen: "<<new_pen<<"   interface pen: "<< ( intP2[1]-intP1[1] ) / ( intP2[0]-intP1[0] ) <<std::endl;
	    }
	  if (CutTypeT2 == 1)
	    {
	      double dxP1 = max_x - x2;
	      double dyP1 = 0 - y2;

	      double dxP2 = tilde_xt2 - x2;
	      double dyP2 = tilde_yt2 - y2;

	      intP1[0] = x2 + sqrt (ratio) * dxP1;
	      intP1[1] = y2 + sqrt (ratio) * dyP1;

	      intP2[0] = x2 + sqrt (ratio) * dxP2;
	      intP2[1] = y2 + sqrt (ratio) * dyP2;
//                 std::cout<<" cut1, start pen: "<<new_pen<<"   interface pen: "<< ( intP2[1]-intP1[1] ) / ( intP2[0]-intP1[0] ) <<std::endl;
	    }
	}

      double point1[2];
      double point2[2];

      if (RotCase == 0)
	{			// point 1 will be at (1,0) with normalized coordinates
	  point1[0] =
	    intP1[0] * cos01 - intP1[1] * sin01 + FullCoordinates[0];
	  point1[1] =
	    intP1[1] * cos01 + intP1[0] * sin01 + FullCoordinates[CNodes + 2];
	  point2[0] =
	    intP2[0] * cos01 - intP2[1] * sin01 + FullCoordinates[0];
	  point2[1] =
	    intP2[1] * cos01 + intP2[0] * sin01 + FullCoordinates[CNodes + 2];
	}
      else
	{			// point 3 will be at (-1,0) with normalized coordinates
	  point1[0] =
	    -intP1[0] * cos03 + intP1[1] * sin03 + FullCoordinates[0];
	  point1[1] =
	    -intP1[1] * cos03 - intP1[0] * sin03 + FullCoordinates[CNodes +
								   2];
	  point2[0] =
	    -intP2[0] * cos03 + intP2[1] * sin03 + FullCoordinates[0];
	  point2[1] =
	    -intP2[1] * cos03 - intP2[0] * sin03 + FullCoordinates[CNodes +
								   2];
	}

      InterfaceCoords.push_back (point1[0]);
      InterfaceCoords.push_back (point2[0]);
      InterfaceCoords.push_back (point1[1]);
      InterfaceCoords.push_back (point2[1]);
    }

  int connect[2];
  int Intcells = InterfaceCoords.size () / (DIMENSION * 2);
  int Intnodes = InterfaceCoords.size () / (DIMENSION);

  if (_ReconInterface != NULL)
    _ReconInterface->decrRef ();
  _ReconInterface = MEDCoupling::MEDCouplingUMesh::New ("RefinedMesh", 1);
  MEDCoupling::DataArrayDouble * CoordArray =
    MEDCoupling::DataArrayDouble::New ();
  CoordArray->alloc (Intnodes, DIMENSION);
  _ReconInterface->allocateCells (Intcells);

  for (int i = 0; i < Intcells; i++)
    {				// REFINED CELLS ================================
      for (int k = 0; k < 2; k++)
	{
	  int new_node_id = k + i * 2;
	  connect[k] = new_node_id;
	  for (int d = 0; d < DIMENSION; d++)
	    {
	      CoordArray->setIJ (connect[k], d,
				 InterfaceCoords[2 * i * DIMENSION + k +
						 d * 2]);
	    }
	}
      _ReconInterface->insertNextCell (INTERP_KERNEL::NORM_SEG2, 2, connect);
    }				// END ADDING REFINED CELLS ==============================================================
  _ReconInterface->finishInsertingCells ();
  _ReconInterface->setCoords (CoordArray);
  _ReconInterface->zipCoords ();

  CoordArray->decrRef ();


  delete[]CellId;
  delete[]dphi;
  delete[]Coordinates;
  delete[]ColVal;
  delete[]FullCoordinates;
  delete[]FullCoordinates_relative;
  return;
}


void
IbUtils::MoveMeshWithStress (MEDCoupling::MEDCouplingFieldDouble *
			     VelocityField, double dt)
{
  std::cout << " -----------------------------------------------\n \
 IBUTILS: updating solid position with viscous stress \n \
----------------------------------------------- \n";

  double Stress[DIMENSION];
  for (int i = 0; i < DIMENSION; i++)
    Stress[i] = 0.;
  ComputeStress (VelocityField, Stress);

  double acceleration[DIMENSION], old_vel[DIMENSION], translation[DIMENSION];

  double gravity[DIMENSION];
  gravity[0] = 0.;
  gravity[1] = -9.81;
#if DIMENSION==3
  gravity[1] = 0;
  gravity[2] = -9.81;
#endif

  for (int i = 0; i < DIMENSION; i++)
    {
      acceleration[i] =
	Stress[i] / (_SolidArea *
		     _rhos) /*+ _buoyancy*gravity[i]* ( 1-_rhof/_rhos ) */ ;
      old_vel[i] = _Velocity[i];
      _Velocity[i] += acceleration[i] * dt;
      translation[i] = 0.5 * (old_vel[i] + _Velocity[i]) * dt;
    }
  _SolidBodyMesh->translate (translation);
  _LinearSolidBodyMesh->translate (translation);

#if DIMENSION==2
  _SolidVelocity->applyFunc (DIMENSION,
			     "IVec*(x+" + to_string (acceleration[0] * dt) +
			     ") + JVec*(y+" +
			     to_string (acceleration[1] * dt) + ")");
#endif
#if DIMENSION==3
  _SolidVelocity->applyFunc (DIMENSION,
			     "IVec*(x+" + to_string (acceleration[0] * dt) +
			     ") + JVec*(y+" +
			     to_string (acceleration[1] * dt) +
			     ")  + KVec*(z+" +
			     to_string (acceleration[2] * dt) + ")");
#endif
  double cd, stress_mod;
  stress_mod = sqrt (Stress[0] * Stress[0] + Stress[1] * Stress[1]);
  cd = stress_mod / (0.5 * _rhof * 0.2 * 0.2 * 0.2);

  UpdateInterpolatedFields (VelocityField);
  return;
}

void
IbUtils::ComputeStress (MEDCoupling::MEDCouplingFieldDouble *
			ReducedVelocityField, double Stress[])
{

  if (_ExtractedVelocity != NULL)
    _ExtractedVelocity->decrRef ();
  _ExtractedVelocity =
    ReducedVelocityField->buildSubPart (_CellToExtract,
					_CellToExtract + _NumCellToExtract);

//     if ( _Proc==0 ) MEDCoupling::WriteField ( "Extracted_vel.med", _ExtractedVelocity, true );
  int Nodes = _ExtractedVelocity->getMesh ()->getNumberOfNodes ();
  int Cells = _ExtractedVelocity->getMesh ()->getNumberOfCells ();
  int Cells2 = _ReconInterface->getNumberOfCells ();
  // velocity components are written as: for i-th node Vel[i*dimension] and Vel[i*dimension+1]
  double *Vel = const_cast < double *>(_ExtractedVelocity->getArray ()->getPointer ());
    
  Stress[0] = 0.;
  Stress[1] = 0.;
  Stress[DIMENSION - 1] = 0.;

  double GlobPressContrib = 0.;
  double GlobViscContrib = 0.;

  double Pa, Pb;
  double Xa[2] = { -0.95, -0.005 };
  double Xb[2] = { -0.85, -0.005 };

  double ViscStress[DIMENSION], PressContrib[DIMENSION];
  for (int dir = 0; dir < DIMENSION; dir++)
    ViscStress[dir] = PressContrib[dir] = 0.;


  const int WeightDim = DIMENSION - 2;
  const int QPoints = _fe[2]->_NoGauss1[WeightDim];

  double InvJac[DIMENSION * DIMENSION], phi_g[3], xyz_g[DIMENSION],
    CanPos[DIMENSION], vel_dx_g[DIMENSION * DIMENSION],
    vel_on_nodes[DIMENSION * NDOF_FEM], Tensor[DIMENSION][DIMENSION],
    Normal[DIMENSION], Tangent[DIMENSION], VelTg[NDOF_FEM],
    Vel_tg_grad_norm_g;

  double dphi[DIMENSION * NDOF_FEM], Pressure[NDOF_P];
  _FamilyType = 1;
  _SrcCoordInterpNodes = 4;

  double NewStress[2];
  NewStress[0] = NewStress[1] = 0.;
  for (int cell = 0; cell < Cells; cell++)
    {				// loop over extracted mesh cells
      std::vector < int >conn;
      double Coord[3 * DIMENSION];
      double MidPoint[DIMENSION];
      MidPoint[0] = 0.;
      MidPoint[1] = 0.;
      MidPoint[DIMENSION - 1] = 0.;
      _ReconInterface->getNodeIdsOfCell (cell, conn);
      for (int dim = 0; dim < conn.size (); dim++)
	{
	  std::vector < double >coord;
	  _ReconInterface->getCoordinatesOfNode (conn[dim], coord);
	  for (int j = 0; j < DIMENSION; j++)
	    {
	      Coord[dim + j * 3] = coord[j];
	      MidPoint[j] += coord[j];
	    }
	}
      Coord[2] = MidPoint[0] / 2;
      Coord[5] = MidPoint[1] / 2;

      std::vector < int >FluidCellConn;
      _ExtractedMesh->getNodeIdsOfCell (cell, FluidCellConn);

      _CoordMatrix.clear ();
      double FEM_coord[DIMENSION * NDOF_FEM];
      for (int s = 0; s < NDOF_FEM; s++)
	{
	  std::vector < double >coord;
	  _ExtractedMesh->getCoordinatesOfNode (FluidCellConn[s], coord);
	  for (int d = 0; d < DIMENSION; d++)
	    {
	      FEM_coord[s + d * NDOF_FEM] = coord[d];
	      if (s < _SrcCoordInterpNodes)
		_CoordMatrix.push_back (coord[d]);
	    }
	}

      double xmin =
	(FEM_coord[0] < FEM_coord[2]) ? FEM_coord[0] : FEM_coord[2];
      double ymin =
	(FEM_coord[0 + NDOF_FEM] <
	 FEM_coord[2 + NDOF_FEM]) ? FEM_coord[0 + NDOF_FEM] : FEM_coord[2 +
									NDOF_FEM];
      double xmax =
	(FEM_coord[0] > FEM_coord[2]) ? FEM_coord[0] : FEM_coord[2];
      double ymax =
	(FEM_coord[0 + NDOF_FEM] >
	 FEM_coord[2 + NDOF_FEM]) ? FEM_coord[0 + NDOF_FEM] : FEM_coord[2 +
									NDOF_FEM];


      for (int dir = 0; dir < DIMENSION; dir++)
	Normal[dir] = _InterfaceNormal[cell * DIMENSION + dir];
      Tangent[0] = Normal[1];
      Tangent[1] = -Normal[0];

      for (int k = 0; k < NDOF_FEM; k++)	// node
	for (int dir = 0; dir < DIMENSION; dir++)	// direction
    vel_on_nodes[k + dir * NDOF_FEM] =
	    Vel[FluidCellConn[k] * (DIMENSION + 1) + dir];
    
      for (int n = 0; n < NDOF_P; n++)
	Pressure[n] = Vel[FluidCellConn[n] * (DIMENSION + 1) + DIMENSION];
    
          for (int nod = 0; nod < NDOF_FEM; nod++)
	{
	  double scal_prod = 0;
	  for (int dim = 0; dim < DIMENSION; dim++)
	    {
	      scal_prod += Tangent[dim] * vel_on_nodes[nod + dim * NDOF_FEM];
	    }
	  int sign = (scal_prod > 0) ? 1 : -1;
	  VelTg[nod] = sign * fabs (scal_prod);
	}

      const double eps = 0.0001;
//                 if((Xa[0]<xmax && Xa[0]>xmin && Xa[1] < ymax && Xa[1] > ymin)){
      _XiEtaChi.clear ();
      XiEtaChiCalc (Xa);
      if (fabs (_XiEtaChi[0]) < 1. + eps && fabs (_XiEtaChi[1]) < 1. + eps)
	{
	  CanPos[0] = _XiEtaChi[0];
	  CanPos[1] = _XiEtaChi[1];
	  Pa = 0.;
	  for (int n = 0; n < NDOF_P; n++)
	    {			// number of dof
	      double lin_phi = LinPhi (n, CanPos);
	      Pa += Pressure[n] * lin_phi;
	    }
	}

//                 if((Xb[0]<xmax && Xb[0]>xmin && Xb[1] < ymax && Xb[1] > ymin)){
      _XiEtaChi.clear ();
      XiEtaChiCalc (Xb);
      if (fabs (_XiEtaChi[0]) < 1. + eps && fabs (_XiEtaChi[1]) < 1. + eps)
	{
	  CanPos[0] = _XiEtaChi[0];
	  CanPos[1] = _XiEtaChi[1];
	  Pb = 0.;
	  for (int n = 0; n < NDOF_P; n++)
	    {			// number of dof
	      double lin_phi = LinPhi (n, CanPos);
	      Pb += Pressure[n] * lin_phi;
	    }
	}


      for (int qp = 0; qp < QPoints; qp++)
	{
	  double det = _fe[2]->JacSur (qp, Coord, InvJac);	// local coord _phi_g and jac
	  double JxW_g = det * _fe[2]->_weight1[WeightDim][qp] /**dt*/ ;	// weight
	  _fe[2]->get_phi_gl_g (DIMENSION - 1, qp, phi_g);	// global coord _phi_g
	  xyz_g[0] = xyz_g[1] = xyz_g[DIMENSION - 1] = 0.;
	  for (int i = 0; i < DIMENSION; i++)
	    for (int j = 0; j < QPoints; j++)
	      xyz_g[i] += Coord[j + i * QPoints] * phi_g[j];

	  _XiEtaChi.clear ();
	  XiEtaChiCalc (xyz_g);
	  CanPos[0] = _XiEtaChi[0];
	  CanPos[1] = _XiEtaChi[1];
	  CanPos[DIMENSION - 1] = _XiEtaChi[DIMENSION - 1];
	  _XiEtaChi.clear ();
	  _fe[2]->get_dphi_on_given_node (DIMENSION, FEM_coord, CanPos, dphi);

	  double p_on_g = 0.;
	  for (int n = 0; n < NDOF_P; n++)
	    {			// number of dof
	      double lin_phi = LinPhi (n, CanPos);
	      p_on_g += Pressure[n] * lin_phi;
	    }
	  // interpolation of velocity derivatives
	  for (int l = 0; l < DIMENSION * DIMENSION; l++)
	    vel_dx_g[l] = 0.;
	  for (int i = 0; i < DIMENSION; i++)	// velocity component  ->  [ux uy vx vy]
	    for (int j = 0; j < DIMENSION; j++)	// direction of derivative
	      for (int n = 0; n < NDOF_FEM; n++)	// number of dof
		vel_dx_g[i * DIMENSION + j] +=
		  vel_on_nodes[n + i * NDOF_FEM] * dphi[n + j * NDOF_FEM];

	  Vel_tg_grad_norm_g = 0;
	  for (int j = 0; j < DIMENSION; j++)	// direction of derivative
	    for (int n = 0; n < NDOF_FEM; n++)	// number of dof
	      Vel_tg_grad_norm_g +=
		VelTg[n] * dphi[n + j * NDOF_FEM] * Normal[j];



	  for (int i = 0; i < DIMENSION; i++)
	    for (int j = 0; j < DIMENSION; j++)
	      Tensor[i][j] =
		/*0.5* */
		(vel_dx_g[i * DIMENSION + j] + vel_dx_g[i + j * DIMENSION]);

	  for (int i = 0; i < DIMENSION; i++)
	    {
	      ViscStress[i] = 0;
	      for (int j = 0; j < DIMENSION; j++)
		ViscStress[i] += JxW_g * _muf * Tensor[i][j] * Normal[j];
	      PressContrib[i] = JxW_g * _rhof * p_on_g * Normal[i];
	    }

	  NewStress[0] +=
	    JxW_g * (Vel_tg_grad_norm_g * Normal[1] * _muf -
		     _rhof * p_on_g * Normal[0]);
	  NewStress[1] +=
	    JxW_g * (Vel_tg_grad_norm_g * Normal[0] * _muf +
		     _rhof * p_on_g * Normal[1]);

	  for (int i = 0; i < DIMENSION; i++)
	    Stress[i] += ViscStress[i] - PressContrib[i];
	  GlobViscContrib += ViscStress[0];
	  GlobPressContrib += PressContrib[0];

	}			// END LOOP OVER GAUSS NODES



    }				// END LOOP OVER CELLS

//     double cd =     NewStress[0] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     double cl = -1.*NewStress[1] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
// //     std::cout<<" ==================================== \n";
// //     std::cout<< "cd:  "<< cd <<"  cl: "<<cl<<std::endl;
// //
//
//     double cd2 =     Stress[0] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     double cl2 = -1.*Stress[1] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     std::cerr<<" ==================================== \n";
//     std::cerr<< "cd:  "<< cd2 <<"  cl: "<<cl2<<std::endl;
//
//     std::cerr<<" ==================================== \n";
  if (_Proc == 0)
    std::cerr << "   " << Pa -
      Pb << "   " << GlobViscContrib << "   " << GlobPressContrib;
  return;
}


void
IbUtils::read_file (int p_femus)
{
  std::string string_value;	// read double, string, dummay
  std::ostringstream file, file1, file2;
  std::vector < std::string > FILES;

  if(p_femus==0) file2 << getenv ("APP_PATH") << "/DATA/MaterialProperties.in";
  else file2 << getenv ("APP_PATH") << "/DATA/DATA"<<p_femus<<"/MaterialProperties.in";

  FILES.push_back (file2.str ());

  for (int i = 0; i < FILES.size (); i++)
    {

      std::ifstream fin;
      fin.open (FILES[i].c_str ());	// stream file
      std::string buf = "";
      if (fin.is_open ())
	{
	  std::cout << "\nInit Reading = " << FILES[i] << std::endl;
	}
      if (fin.is_open ())
	{
	  while (buf != "/")
	    {
	      fin >> buf;	// find "/" file start
	    }
	  fin >> buf;
	  while (buf != "/")
	    {
	      if (buf == "#")
		{
		  getline (fin, buf);	// comment line
		}
	      else
		{
		  fin >> string_value;
		  _FileMap[buf] = string_value;
		  std::cout << buf << "\t" << string_value << std::endl;
		}
	      fin >> buf;
	    }
	}
      else
	{
	  std::cerr << "IBMover.read_file(): no parameter file found\t ->" <<
	    FILES[i] << std::endl;
	  abort ();
	}
      std::cout << "IBMover.read_file() End Reading file " << FILES[i] <<
	std::endl;
      fin.close ();
    }


  if (_FileMap["rho0"] != "")
    {
      _rhof = stof (_FileMap["rho0"]);
    }
  else
    {
      _rhof = 1;
      std::cout <<
	" IBMover: Fluid density not found. Assuming rho_fluid = 1 \n";
      std::cerr <<
	" IBMover: Fluid density not found. Assuming rho_fluid = 1 \n";
    }
  if (_FileMap["mu0"] != "")
    {
      _muf = stof (_FileMap["mu0"]);
    }
  else
    {
      _muf = 1;
      std::cout <<
	" IBMover: Fluid viscosity not found. Assuming mu_fluid = 1 \n";
      std::cerr <<
	" IBMover: Fluid viscosity not found. Assuming mu_fluid = 1 \n";
    }
  if (_FileMap["buoyancy"] != "")
    {
      _buoyancy = stof (_FileMap["buoyancy"]);
    }
  else
    {
      _buoyancy = 0;
      std::cout <<
	" IBMover: Buoyancy flag not found. Assuming no buoyancy \n";
      std::cerr <<
	" IBMover: Buoyancy flag not found. Assuming no buoyancy \n";
    }
  if (_FileMap["rhosolid"] != "")
    {
      _rhos = stof (_FileMap["rhosolid"]);
    }
  else
    {
      _rhos = _rhof;
      std::cout <<
	" IBMover: Solid density not found. Assuming rho_solid equal to rho_fluid \n";
      std::cerr <<
	" IBMover: Solid density not found. Assuming rho_solid equal to rho_fluid \n";
    }
  if (_FileMap["LowerThreshold"] != "")
    {
      _LowerThreshold = stof (_FileMap["LowerThreshold"]);
    }
  else
    {
      _LowerThreshold = 0.001;
      std::cout <<
	" IBMover: LowerThreshold not found. Assuming default 0.001 \n";
      std::cerr <<
	" IBMover: LowerThreshold not found. Assuming default 0.001 \n";
    }
  if (_FileMap["UpperThreshold"] != "")
    {
      _UpperThreshold = stof (_FileMap["UpperThreshold"]);
    }
  else
    {
      _UpperThreshold = 0.999;
      std::cout <<
	" IBMover: UpperThreshold not found. Assuming default 0.999 \n";
      std::cerr <<
	" IBMover: UpperThreshold not found. Assuming default 0.999 \n";
    }
  if (_FileMap["SimpleIBMethod"] != "")
    {
      _SimpleIBMethod = stoi (_FileMap["SimpleIBMethod"]);
    }
  else
    {
      _SimpleIBMethod = 0;
      std::cout <<
	" IBMover: SimpleIBMethod not found. Assuming default 0 \n";
      std::cerr <<
	" IBMover: SimpleIBMethod not found. Assuming default 0 \n";
    }
  return;

}

void
IbUtils::InitSolidVelocity (vector < double >Velocity)
{
  int dim = Velocity.size ();
  for (int i = 0; i < dim; i++)
    _Velocity[i] = Velocity[i];
  std::string VelField =
    "IVec*" + std::to_string (Velocity[0]) + " + JVec*" +
    std::to_string (Velocity[1]);
  if (dim == 3)
    VelField += " + KVec*" + std::to_string (Velocity[2]);
  _SolidVelocity->applyFunc (dim, VelField);
  return;
}

void
IbUtils::InitPosition (double FirstTranslation[],
		       MEDCoupling::MEDCouplingFieldDouble * VelocityField)
{
  _SolidBodyMesh->translate (FirstTranslation);
  _LinearSolidBodyMesh->translate (FirstTranslation);
  UpdateInterpolatedFields (VelocityField);
  return;
}

void
IbUtils::MoveWithGivenVectorAndRotation (double AngVel,
					 double center[],
					 double vector[],
					 double dt,
					 double velocity[],
					 MEDCoupling::MEDCouplingFieldDouble *
					 VelocityField)
{
  double Translation[3] =
    { velocity[0] * dt, velocity[1] * dt, velocity[2] * dt };
  _SolidBodyMesh->translate (Translation);
  _SolidBodyMesh->rotate (center, vector, AngVel * dt);
  _LinearSolidBodyMesh->translate (Translation);
  _LinearSolidBodyMesh->rotate (center, vector, AngVel * dt);

  std::vector < double >Coordinates;
  double dist, mod, norm_x, norm_y, v_x, v_y;
  const int Nnodes = _SolidBodyMesh->getNumberOfNodes ();
  for (int Node = 0; Node < Nnodes; Node++)
    {
      _SolidBodyMesh->getCoordinatesOfNode (Node, Coordinates);
      dist =
	sqrt (Coordinates[0] * Coordinates[0] +
	      Coordinates[1] * Coordinates[1] + 1.e-20);
      mod = dist * (AngVel);
      norm_x = ((-1) * Coordinates[1]) / dist;
      norm_y = Coordinates[0] / dist;
      v_x = norm_x * mod + velocity[0];
      v_y = norm_y * mod + velocity[1];
      _SolidVelocity->getArray ()->setIJ (Node, 0, v_x);
      _SolidVelocity->getArray ()->setIJ (Node, 1, v_y);
      Coordinates.clear ();
    }

  UpdateInterpolatedFields (VelocityField);
  return;
}

void
IbUtils::UpdateInterpolatedFields (MEDCoupling::MEDCouplingFieldDouble *
				   VelocityField)
{

  if (_FluidCellColor != NULL)
    _FluidCellColor->decrRef ();
  if (_FluidNodeColor != NULL)
    _FluidNodeColor->decrRef ();
  if (_InterpolatedSolidVelocity != NULL)
    _InterpolatedSolidVelocity->decrRef ();

  // PIECEWISE VOLUME FRACTION
  MEDCoupling::MEDCouplingRemapper Remap;
  Remap.setPrecision (1e-11);
  Remap.setIntersectionType (_PieceWiseInterpScheme);	// <- more accurate than triangulate
  Remap.prepare (_LinearSolidBodyMesh, _LinearFluidMesh, "P0P0");
  _FluidCellColor = Remap.transferField (_SolidCellColor, 0.);

  _RemapMatrix.clear ();
  _RemapMatrix = Remap.getCrudeMatrix ();
  ofstream CoolerTemperature;
  CoolerTemperature.open ("CrudeMatrix.dat");
  CoolerTemperature <<
    "Target Mesh Cell # n : (Source Mesh Cell ID, Overlapping Area)\n====================================================\n";

  int cells = _FluidCellColor->getMesh ()->getNumberOfCells ();
  for (int i = 0; i < cells; i++)
    {
      double col = _FluidCellColor->getArray ()->getIJ (i, 0);
      if (col > 1.1)
	{
	  std::cout << "\033[1;21m IBMOVER: Col > 1 " << col << "  cell " << i
	    << "\033[0m\n";
	for (auto & x:_RemapMatrix[i])
	    CoolerTemperature << "(" << x.first << ", " << x.second << "),\t";
	  CoolerTemperature << "\n";
	  _FluidCellColor->getArray ()->setIJ (i, 0, 0.5 * col);
	}
    }

  CoolerTemperature.close ();
  // NODAL VOLUME FRACTION AND VELOCITY FIELD
  std::cout << "\033[1;21m IBMOVER: FILL PARAMETERS\033[0m\n";
  FillParameters (_SolidBodyMesh, _FluidMesh, Volume, 1.e-5);
  _FluidNodeColor = InterpolatedField (_SolidNodeColor, 0., 2);
  _InterpolatedSolidVelocity = InterpolatedField (_SolidVelocity, VelocityField);

  return;
}

void
IbUtils::UpdateFluidMesh (const MEDCoupling::MEDCouplingUMesh * NewMesh)
{

  _FluidMesh->decrRef ();
  _LinearFluidMesh->decrRef ();

  _FluidMesh = NewMesh->deepCopy ();

  MEDCoupling::MEDCouplingUMesh * FemusMesh2 = NewMesh->deepCopy ();
  FemusMesh2->convertQuadraticCellsToLinear ();
  _LinearFluidMesh = FemusMesh2->deepCopy ();

  FemusMesh2->decrRef ();

  FillParameters (_SolidBodyMesh, _FluidMesh, Volume, 1.e-5);
  
  return;
}


void
IbUtils::UpdateSolidMesh (const MEDCoupling::MEDCouplingUMesh * NewMesh)
{

  _SolidBodyMesh->decrRef ();
  _LinearSolidBodyMesh->decrRef ();

  _SolidBodyMesh = NewMesh->deepCopy ();

  MEDCoupling::MEDCouplingUMesh * FemusMesh2 = NewMesh->deepCopy ();
  FemusMesh2->convertQuadraticCellsToLinear ();
  _LinearSolidBodyMesh = FemusMesh2->deepCopy ();

  FemusMesh2->decrRef ();

  FillParameters (_FluidMesh, _SolidBodyMesh, Volume, 1.e-5);
  
  return;
}


void IbUtils::UpdateSolidVelocity (const MEDCoupling::MEDCouplingFieldDouble * SourceField)
{
     
  _SolidVelocity->decrRef ();
  _SolidVelocity = SourceField->deepCopy ();
  
  return;
    
}


void
IbUtils::InsertTriangulatedCell (int SourceCellId,
				 std::vector < int >&NodesPerCell,
				 std::vector < int >OrderedNodes,
				 std::vector < int >&Sequence,
				 std::vector < double >&TriangulationCoords,
				 double FullCoordinates[], int offset)
{

  if (Sequence.size () == 3)
    {

      double vec_prod =
	(FullCoordinates[OrderedNodes[Sequence[1]]] -
	 FullCoordinates[OrderedNodes[Sequence[0]]]) *
	(FullCoordinates[OrderedNodes[Sequence[2]] + offset] -
	 FullCoordinates[OrderedNodes[Sequence[1]] + offset]) -
	(FullCoordinates[OrderedNodes[Sequence[1]] + offset] -
	 FullCoordinates[OrderedNodes[Sequence[0]] +
			 offset]) *
	(FullCoordinates[OrderedNodes[Sequence[2]]] -
	 FullCoordinates[OrderedNodes[Sequence[1]]]);

      if (vec_prod < 0)
	{
	  double a = Sequence[1], b = Sequence[2];
	  Sequence[1] = b;
	  Sequence[2] = a;
	}
    }



  _TriangulatedCellsPerSourceCell[SourceCellId].
    push_back (NodesPerCell.size ());
  NodesPerCell.push_back (NodesPerCell[NodesPerCell.size () - 1] +
			  Sequence.size ());
  for (int l = 0; l < Sequence.size (); l++)
    for (int dim = 0; dim < DIMENSION; dim++)
      TriangulationCoords.push_back (FullCoordinates
				     [OrderedNodes[Sequence[l]] +
				      dim * offset]);

  return;
}


void
IbUtils::GetCellContainingPoint (int TargetCellId, double NodeCoord[],
				 double CanPos[], int &CellId)
{

  std::vector < int >SourceIDs;
for (auto & x:_RemapMatrix[TargetCellId])
    SourceIDs.push_back (x.first);

  int nCells = SourceIDs.size ();
  int *IDs = new int[nCells];
  for (int i = 0; i < nCells; i++)
    IDs[i] = SourceIDs[i];

  MEDCoupling::MEDCouplingFieldDouble * Part =
    _SolidVelocity->buildSubPart (IDs, IDs + nCells);


  MEDCoupling::MCAuto < MEDCoupling::DataArrayInt > elts;
  MEDCoupling::MCAuto < MEDCoupling::DataArrayInt > eltsIndex;
  Part->getMesh ()->getCellsContainingPoints (NodeCoord, 1, 1.e-5, elts,
					      eltsIndex);

  int NumPossibleCells = eltsIndex->getIJ (1, 0) - eltsIndex->getIJ (0, 0);
  bool found = false;

  std::vector < int >SourceConn;
  if (NumPossibleCells > 0)
    {
      int iCount = 0;
      while (!found && iCount < NumPossibleCells)
	{

	  _XiEtaChi.clear ();
	  // cell ids are elts[ eltsIndex[ i ]],..., elts[ eltsIndex[ i ] + NumPossibleCells ].
	  int NewCell = elts->getIJ (eltsIndex->getIJ (0, 0) + iCount, 0);
	  CellId = NewCell;
	  Part->getMesh ()->getNodeIdsOfCell (NewCell, SourceConn);

	  for (int SNode = 0; SNode < _SrcCellNodes; SNode++)
	    {
	      if (SNode < _SrcCoordInterpNodes)
		{
		  Part->getMesh ()->getCoordinatesOfNode (SourceConn[SNode], _CoordMatrix);	// only linear nodes
		}
	    }
	  XiEtaChiCalc (NodeCoord);
	  int Contained = IsNodeContained ();
	  if (Contained == 0)
	    {
	      found = false;
	      _XiEtaChi.clear ();
	      _CoordMatrix.clear ();
	      SourceConn.clear ();
	    }
	  else if (Contained == 1)
	    found = true;
	  iCount++;
	}			// end while

      if (found)
	{
//                     CellId = NewCell;
	  for (int Mcomp = 0; Mcomp < _MeshDim; Mcomp++)
	    {
	      CanPos[Mcomp] = _XiEtaChi[Mcomp];
	    }
	}
      _CoordMatrix.clear ();
      _XiEtaChi.clear ();
    }				// end  if cell==-1
  if (!found)
    CellId = -1;
  SourceConn.clear ();

  Part->decrRef ();

  return;
}

double
IbUtils::TriangleArea (double *point1, double *point2, double *point3)
{
  double Area = 0;

  double det =
    point1[0] * point2[1] + point2[0] * point3[1] + point3[0] * point1[1] -
    (point1[0] * point3[1] + point3[0] * point2[1] + point2[0] * point1[1]);

  Area = 0.5 * fabs (det);
  return Area;
}

double
IbUtils::TriangleArea (double *points /* x1 x2 x3, y1 y2 y3 */ )
{
  double Area = 0.;
  double det = 0.;
  for (int i = 0; i < 3; i++)
    {
      int col_plus = (i + 1) % 3;
      int col_minus = (i + 2) % 3;
      det += points[i] * (points[3 + col_plus] - points[3 + col_minus]);
    }
  Area = 0.5 * fabs (det);
  return Area;
}

double
IbUtils::QuadrangleArea (double *points /* x1 x2 x3 x4, y1 y2 y3 y4 */ )
{
  double Area = 0.;
  double det = 0.;
  for (int i = 0; i < 4; i++)
    {
      int col_plus = (i + 1) % 4;
      int col_minus = (i + 3) % 4;
      det += points[i] * (points[4 + col_plus] - points[4 + col_minus]);
    }
  Area = 0.5 * fabs (det);
  return Area;
}


MEDCoupling::MEDCouplingFieldDouble * IbUtils::InterpolateSolidOnFluid (MEDCoupling::MEDCouplingFieldDouble *
				   SolidVelocity, MEDCoupling::MEDCouplingFieldDouble * FluidVelocity)
{
  MEDCoupling::MEDCouplingFieldDouble * f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  // NODAL VOLUME FRACTION AND VELOCITY FIELD
  std::cout << "\033[1;21m IBMOVER: FILL PARAMETERS\033[0m\n";
  FillParameters (_SolidBodyMesh, _FluidMesh, Volume, 1.e-5);
  f = InterpolatedField (SolidVelocity, FluidVelocity);
  
  return f;
}


MEDCoupling::MEDCouplingFieldDouble * IbUtils::InterpolateFluidOnSolid (MEDCoupling::MEDCouplingFieldDouble *
				   VelocityField, MEDCoupling::MEDCouplingFieldDouble * SolidBoundaryField, 
                   const MEDCoupling::MEDCouplingUMesh *BoundaryMesh)
{
  MEDCoupling::MEDCouplingFieldDouble * f = MEDCoupling::MEDCouplingFieldDouble::New ( MEDCoupling::ON_NODES );
  // NODAL VOLUME FRACTION AND VELOCITY FIELD
  std::cout << "\033[1;21m IBMOVER: FILL PARAMETERS\033[0m\n";
  FillParameters (_FluidMesh, BoundaryMesh, Volume, 1.e-5);
  f = InterpolatedField (VelocityField, SolidBoundaryField);
  
  return f;
}


MEDCoupling::MEDCouplingFieldDouble * IbUtils::ComputeStressField (MEDCoupling::MEDCouplingFieldDouble *
			ReducedVelocityField, const MEDCoupling::MEDCouplingUMesh *BoundaryMesh )
{

//   if (_ExtractedVelocity != NULL)
//     _ExtractedVelocity->decrRef ();
//   _ExtractedVelocity =
//     ReducedVelocityField->buildSubPart (_CellToExtract,
// 					_CellToExtract + _NumCellToExtract);

// //     if ( _Proc==0 ) MEDCoupling::WriteField ( "Extracted_vel.med", _ExtractedVelocity, true );
  int Nodes = ReducedVelocityField->getMesh ()->getNumberOfNodes ();
  int Cells = ReducedVelocityField->getMesh ()->getNumberOfCells ();
  // velocity components are written as: for i-th node Vel[i*dimension] and Vel[i*dimension+1]
  double *Vel = const_cast < double *>(ReducedVelocityField->getArray ()->getPointer ());
  double Stress[DIMENSION];
  std::vector<double> StressVector( Nodes * DIMENSION );
//   StressVector.resize(Nodes * DIMENSION);
    
  Stress[0] = 0.;
  Stress[1] = 0.;
  Stress[DIMENSION - 1] = 0.;

  double GlobPressContrib = 0.;
  double GlobViscContrib = 0.;

  double Pa, Pb;
  double Xa[2] = { -0.95, -0.005 };
  double Xb[2] = { -0.85, -0.005 };

  double ViscStress[DIMENSION], PressContrib[DIMENSION];
  for (int dir = 0; dir < DIMENSION; dir++)
    ViscStress[dir] = PressContrib[dir] = 0.;


  const int WeightDim = DIMENSION - 2;
  const int QPoints = _fe[2]->_NoGauss1[WeightDim];

  double InvJac[DIMENSION * DIMENSION], phi_g[3], xyz_g[DIMENSION],
    CanPos[DIMENSION], vel_dx_g[DIMENSION * DIMENSION],
    vel_on_nodes[DIMENSION * NDOF_FEMB], Tensor[DIMENSION][DIMENSION],
    Normal[DIMENSION], Tangent[DIMENSION], VelTg[NDOF_FEMB],
    Vel_tg_grad_norm_g;

  double dphi[DIMENSION * NDOF_FEMB], Pressure[NDOF_PB];
  _FamilyType = 1;
  _SrcCoordInterpNodes = 4;

  double NewStress[2];
  NewStress[0] = NewStress[1] = 0.;
  
    MEDCoupling::MEDCouplingFieldDouble * StressFieldBD;
    StressFieldBD = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES);
    StressFieldBD->setMesh(BoundaryMesh);    
    StressFieldBD->setName("Stress");    
    StressFieldBD->setNature(MEDCoupling::IntensiveMaximum);
    
    MEDCoupling::DataArrayDouble *array = MEDCoupling::DataArrayDouble::New();
    array->alloc(Nodes,DIMENSION);    array->fillWithZero();

  for (int cell = 0; cell < Cells; cell++)
    {				// loop over extracted mesh cells
      std::vector < int >conn;
      double Coord[3 * DIMENSION];
      BoundaryMesh->getNodeIdsOfCell (cell, conn);
      for (int dim = 0; dim < conn.size (); dim++)
	{
	  std::vector < double >coord;
	  BoundaryMesh->getCoordinatesOfNode (conn[dim], coord);
	  for (int j = 0; j < DIMENSION; j++)
	    {
	      Coord[dim + j * 3] = coord[j];
	    }
	}

      std::vector < int >FluidCellConn;
//       _ExtractedVelocity->getMesh ()->getNodeIdsOfCell (cell, FluidCellConn);
      BoundaryMesh->getNodeIdsOfCell (cell, FluidCellConn);

//     for (int dir = 0; dir < DIMENSION; dir++)   Normal[dir] = _InterfaceNormal[cell * DIMENSION + dir];
    _fe[2]->normal_g (Coord, Normal);
    Tangent[0] = Normal[1];
    Tangent[1] = -Normal[0];

    for (int k = 0; k < NDOF_FEMB; k++)	// node
        for (int dir = 0; dir < DIMENSION; dir++)	// direction
            vel_on_nodes[k + dir * NDOF_FEMB] =
            Vel[FluidCellConn[k] * (DIMENSION + 1) + dir];
    
    for (int n = 0; n < NDOF_PB; n++)
        Pressure[n] = Vel[FluidCellConn[n] * (DIMENSION + 1) + DIMENSION];
    
    for (int nod = 0; nod < NDOF_FEMB; nod++)
    {
        double scal_prod = 0;
        for (int dim = 0; dim < DIMENSION; dim++)
	    {
            scal_prod += Tangent[dim] * vel_on_nodes[nod + dim * NDOF_FEMB];
	    }
        int sign = (scal_prod > 0) ? 1 : -1;
        VelTg[nod] = sign * fabs (scal_prod);
	}

      const double eps = 0.0001;
//                 if((Xa[0]<xmax && Xa[0]>xmin && Xa[1] < ymax && Xa[1] > ymin)){
    _XiEtaChi.clear ();
    XiEtaChiCalc (Xa);
    if (fabs (_XiEtaChi[0]) < 1. + eps && fabs (_XiEtaChi[1]) < 1. + eps)
    {
        CanPos[0] = _XiEtaChi[0];
        CanPos[1] = _XiEtaChi[1];
        Pa = 0.;
        for (int n = 0; n < NDOF_PB; n++)
	    {			// number of dof
            double lin_phi = LinPhi (n, CanPos);
            Pa += Pressure[n] * lin_phi;
	    }
	}

//                 if((Xb[0]<xmax && Xb[0]>xmin && Xb[1] < ymax && Xb[1] > ymin)){
      _XiEtaChi.clear ();
      XiEtaChiCalc (Xb);
    if (fabs (_XiEtaChi[0]) < 1. + eps && fabs (_XiEtaChi[1]) < 1. + eps)
    {
        CanPos[0] = _XiEtaChi[0];
        CanPos[1] = _XiEtaChi[1];
        Pb = 0.;
        for (int n = 0; n < NDOF_PB; n++)
	    {			// number of dof
            double lin_phi = LinPhi (n, CanPos);
            Pb += Pressure[n] * lin_phi;
	    }
	}

    for (int qp = 0; qp < QPoints; qp++)
    {
        double det = _fe[2]->JacSur (qp, Coord, InvJac);	// local coord _phi_g and jac
        double JxW_g = det * _fe[2]->_weight1[WeightDim][qp] /**dt*/ ;	// weight
        _fe[2]->get_phi_gl_g (DIMENSION - 1, qp, phi_g);	// global coord _phi_g
        xyz_g[0] = xyz_g[1] = xyz_g[DIMENSION - 1] = 0.;
        for (int i = 0; i < DIMENSION; i++)
            for (int j = 0; j < QPoints; j++)
                xyz_g[i] += Coord[j + i * QPoints] * phi_g[j];

	  _XiEtaChi.clear ();
	  XiEtaChiCalc (xyz_g);
	  CanPos[0] = _XiEtaChi[0];
	  CanPos[1] = _XiEtaChi[1];
	  CanPos[DIMENSION - 1] = _XiEtaChi[DIMENSION - 1];
	  _XiEtaChi.clear ();
      _fe[2]->get_dphi_gl_g ( DIMENSION-1,qp,InvJac,dphi );       // global coord deriv
	  double p_on_g = 0.;
	  for (int n = 0; n < NDOF_PB; n++)
	    {			// number of dof
	      double lin_phi = LinPhi (n, CanPos);
	      p_on_g += Pressure[n] * lin_phi;
	    }
	  // interpolation of velocity derivatives
	  for (int l = 0; l < DIMENSION * DIMENSION; l++)
	    vel_dx_g[l] = 0.;
	  for (int i = 0; i < DIMENSION; i++)	// velocity component  ->  [ux uy vx vy]
	    for (int j = 0; j < DIMENSION; j++)	// direction of derivative
	      for (int n = 0; n < NDOF_FEMB; n++)	// number of dof
		vel_dx_g[i * DIMENSION + j] +=
		  vel_on_nodes[n + i * NDOF_FEMB] * dphi[n + j * NDOF_FEMB];

	  Vel_tg_grad_norm_g = 0;
	  for (int j = 0; j < DIMENSION; j++)	// direction of derivative
	    for (int n = 0; n < NDOF_FEMB; n++)	// number of dof
	      Vel_tg_grad_norm_g +=
		VelTg[n] * dphi[n + j * NDOF_FEMB] * Normal[j];



	  for (int i = 0; i < DIMENSION; i++)
	    for (int j = 0; j < DIMENSION; j++)
	      Tensor[i][j] =
		/*0.5* */
		(vel_dx_g[i * DIMENSION + j] + vel_dx_g[i + j * DIMENSION]);

	  for (int i = 0; i < DIMENSION; i++)
	    {
	      ViscStress[i] = 0;
	      for (int j = 0; j < DIMENSION; j++)  ViscStress[i] += JxW_g * _muf * Tensor[i][j] * Normal[j];
	      PressContrib[i] = JxW_g * _rhof * p_on_g * Normal[i];
	    }

	  NewStress[0] +=
	    JxW_g * (Vel_tg_grad_norm_g * Normal[1] * _muf -
		     _rhof * p_on_g * Normal[0]);
	  NewStress[1] +=
	    JxW_g * (Vel_tg_grad_norm_g * Normal[0] * _muf +
		     _rhof * p_on_g * Normal[1]);
	  for (int i = 0; i < DIMENSION; i++){
	    Stress[i] = ViscStress[i] - PressContrib[i];
        StressVector[FluidCellConn[qp] * (DIMENSION) + i] = Stress[i];
        array->setIJ(FluidCellConn[qp],i,Stress[i]);
      }
	  GlobViscContrib += ViscStress[0];
	  GlobPressContrib += PressContrib[0];
      

	}			// END LOOP OVER GAUSS NODES


    }				// END LOOP OVER CELLS
    

//     double cd =     NewStress[0] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     double cl = -1.*NewStress[1] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
// //     std::cout<<" ==================================== \n";
// //     std::cout<< "cd:  "<< cd <<"  cl: "<<cl<<std::endl;
// //
//
//     double cd2 =     Stress[0] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     double cl2 = -1.*Stress[1] * 2 /(_rhof * 0.1 * 0.2 * 0.2);
//     std::cerr<<" ==================================== \n";
//     std::cerr<< "cd:  "<< cd2 <<"  cl: "<<cl2<<std::endl;
//
//     std::cerr<<" ==================================== \n";
    if (_Proc == 0)
        std::cerr << "   " << Pa - Pb << "   " << GlobViscContrib << "   " << GlobPressContrib;
     
    for(int kj=0; kj<DIMENSION; kj++)   array->setInfoOnComponent(kj,std::to_string(kj)); // set info -> array
    StressFieldBD->setArray(array);    // set array -> f
    array->decrRef();     // delete array
    StressFieldBD->checkConsistencyLight();  // check f
        
    StressVector.clear();
      
  return StressFieldBD;
}

MEDCoupling::MEDCouplingFieldDouble * IbUtils::GetFluidCellColor ()
{
  return _FluidCellColor->deepCopy ();
}

MEDCoupling::MEDCouplingFieldDouble * IbUtils::GetFluidNodeColor ()
{
  return _FluidNodeColor->deepCopy ();
}

MEDCoupling::MEDCouplingFieldDouble * IbUtils::GetInterpolatedSolidVelocity ()
{
  return _InterpolatedSolidVelocity->deepCopy ();
}
