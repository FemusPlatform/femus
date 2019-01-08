#ifndef __BoundInterp__
#define __BoundInterp__

#include <vector>
#include "MMed.h"
#include "Solverlib_conf.h"
#include<map>
#ifdef HAVE_MED
namespace MEDCoupling {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
class DataArrayInt;
class DataArrayDouble;
}

//! Method for the integration function #getMediumValuesOnBoundary_elem
// enum method{
//    Mean=0,    /**< Mean integral value */
//    AxiMean=1, /**< Mean integral value with axisymmetry */
//    Bulk=2,    /**< Mean integral value weighted on velocity component normal to surface */
//    AxiBulk=3, /**< Mean integral value weighted on velocity component normal to surface and axisymmetry */
//    Area=4,    /**< Area of the integration surface */
//    Integral=5 /**< Integral value */
// };  

//! Domain type of mesh group
/*!
 * Domain type of the mesh where the solution interpolation is performed.
 * We need to distinguish the interpolation on boundary groups and on volume groups because
 * the algorithm for searching the target node inside the source mesh is different
 */
//  enum DomainType{Boundary=0, Volume=1};

/// Class for the interpolation of a solution from a source mesh to a target mesh 
class BoundInterp : public MMed {  
public:

   BoundInterp ();
   
   //! Constructor of the BoundInterp Class
   BoundInterp  (
           const MEDCoupling::MEDCouplingUMesh * SourceMesh, /**< Mesh support of the source geometry */
	       const MEDCoupling::MEDCouplingUMesh * TargetMesh, /**< Mesh support of the target geometry */
	       DomainType bdd= Boundary /**< Domain type of the mesh group (Boundary of Volume) */
	      );
   BoundInterp  (
           const MEDCoupling::MEDCouplingUMesh * SourceMesh, /**< Mesh support of the source geometry */
	       const MEDCoupling::MEDCouplingUMesh * TargetMesh, /**< Mesh support of the target geometry */
           int procId,
            DomainType bdd = Boundary /**< Domain type of the mesh group (Boundary of Volume) */
	      );

   //! Destructor of the BoundInterp Class
  ~BoundInterp();
   void terminate();
   void setMesh();
   
   void FillParameters(
     const MEDCoupling::MEDCouplingUMesh * SourceMesh, ///< Mesh support of the source geometry 
     const MEDCoupling::MEDCouplingUMesh * TargetMesh, ///< Mesh support of the target geometry 
     int DomainType = Boundary,                        ///< Domain type of the mesh group (Boundary of Volume) 
     double XiEtaToll = 1.e-3                          ///< tolerance
  );
   
   //! Space Dimension of the solved geometry
   int _SpaceDim;
   //! triangular(0) quadrangular(1)
   int _FamilyType;
   
   //! Mesh Dimension
   /*!
    * Dimension of the mesh group: if the group is a volume-one then _MeshDim = _SpaceDim,
    * while if the group is a boundary-one then _MeshDim = _SpaceDim -1
    */
   int _MeshDim;
   
   //! Number of nodes per cell in source mesh
   /*! Number of nodes (per cell) used in the source mesh for the interpolation of the coordinates
    */
   int _SrcCoordInterpNodes;
   
   //! Number of nodes per cell in source mesh
   /*! Number of nodes (per cell) used in the source mesh for the interpolation of the solution
    * This value is different to #_SrcCoordInterpNodes because _SrcCoordInterpNodes are for linear interpolation
    * while _SrcCellNodes are for quadratic interpolation
    */
   int _SrcCellNodes;
   
   //! Total number of nodes in the target mesh
   int _TrgNodes;
   
   //! Total number of cells in the target mesh
   int _TrgCells;
   
   //! Total number of cells in the source mesh
   int _SrcCells;
   int _AlreadyInitialized;
   
   //! Vector where the gradient components of function F are stored
   std::vector<double> _dF;
   
   //! Vector containing the Hessian matrix components for function F
   std::vector<double> _d2F;
   
   //! Vector containing coordinates of source mesh vertex nodes {xyz xyz xyz}
   std::vector<double> _CoordMatrix; 
   
   //! Vector of \f$ (\xi,\eta,\chi) \f$
   std::vector<double> _XiEtaChi;
   //! Check
   int IsNodeContained ();
   
   
   
   // ====================================================================================================================
   //! Nodes of source mesh
   /*!
    * This is a MEDCoupling::DataArrayInt* where we store the ids of the source mesh nodes 
    * that are used for the interpolation of the solution on a given target mesh node.
    * This object contains #_TrgNodes tuples, each one with dimension #_SrcCellNodes
    */
   MEDCoupling::DataArrayInt * _BoundingNodes=NULL;
   // ====================================================================================================================
   //! Coordinates in reference element
   /*!
    * This is a MEDCoupling::DataArrayDouble* where we store the coordinates in the reference
    * element of a given point of the target mesh.that are used for the interpolation of the solution on a given target mesh node.
    * This object contains #_TrgNodes tuples, each one with dimension #_MeshDim
    */
   MEDCoupling::DataArrayDouble * _XiEta=NULL;
   // ====================================================================================================================
   //! \f$ (\xi,\eta) \f$ coordinates for boundary mesh groups
   /*!
    * The function returns the coordinates \f$ (\xi,\eta) \f$ (if _SpaceDim=3) 
    * or \f$ (\xi) \f$ (if _SpaceDim=2) for boundary mesh groups.
    * The coordinates ar from the XiEta DataArrayDouble* and are written
    * inside the XiEta vector given as input parameter for the GetXiEta function
    */
   void GetXiEta(int node, /**< Id of target mesh node */
		 std::vector<double> &XiEta /**< Vector where the coordinates of source mesh node in reference element are stored */
		);
   // ====================================================================================================================
   //! \f$ (\xi,\eta,\chi) \f$ coordinates for volume mesh groups
   /*!
    * The function returns the coordinates \f$ (\xi,\eta,\chi) \f$ (if _SpaceDim=3) 
    * or \f$ (\xi,\eta) \f$ (if _SpaceDim=2) for volume mesh groups.
    * The coordinates are read from the XiEta DataArrayDouble* and are written
    * inside the XiEtaChi vector given as input parameter for the GetXiEtaChi function
    */
   void GetXiEtaChi(int node, /**< Id of target mesh node */
		    std::vector<double> &XiEtaChi /**< Vector where the coordinates of source mesh node in reference element are stored */
		   );
   // ====================================================================================================================
   //! Source mesh interpolating nodes
   /*!
    * The function returns the source mesh nodes that are used for interpolationg the solution on target mesh node <node>.
    * The nodes are read from the #BoundingNodes DataArrayInt* and written into the VectorContainingNodes vector given
    * as input parameter for the GetInterNodes function
    */
   void GetInterNodes(int node, /**< Id of target mesh node */
		      std::vector<int> &VectorContainingNodes /**< Vector where the ids of source mesh nodes are stored */
		     );
   

// ====================================================================================================================
 //! Function for calculation of \f$ (\xi,\eta,\chi) \f$ for a surface
   /*!
    * The function is called when the calcualtion of reference element coordinates 
    * is performed on a boundary-type mesh group.
    * The coordinate values are stored inside the XiEtaBound vector
    */
void XiEtaCalc_2D(
  double NodePos[],                // physical coord (in)
  std::vector<double> &XiEtaBound,  // reference coord (out)
   const int  Quad [],
   int npt_el=4
);
// ================================================================================================
 //! Function for calculation of \f$ (\xi,\eta,\chi) \f$ for a line
   /*!
    * The function is called when the calcualtion of reference element coordinates 
    * is performed on a boundary-type mesh group.
    * The coordinate values are stored inside the XiEtaBound vector
    */
void XiEtaCalc_1D(
  double NodePos[],                // physical coord (in)
  std::vector<double> &XiEtaBound,  // reference coord (out)
int  npt_el=2
);
// ================================================================================================
   //! Function for calculation of \f$ (\xi,\eta,\chi) \f$
   /*!
    * The function is called when the calcualtion of reference element coordinates is performed on a volume-type mesh group.
    * The coordinate values are stored inside the #_XiEtaChi vector
    */
   void XiEtaChiCalc(
      double NodePos[] /**< Coordinates of the target mesh node */
   );
   
  // ================================================================================================ 
   //! Squared error distance for \f$ (\xi,\eta,\chi) \f$ calculation.
   /*!
     This function is used when the user is interested into a volume-type interpolation of the solution (surface in 2D and volume in 3D).
     The calculation of the node \f$ \vec{x}_P \f$ coordinates in the reference element is performed iteratively. 
     Given the tuple \f$ (\xi,\eta,\chi) \f$ it computes the quantity 
     \f$ \\ F = \sum\limits_{i=1}^{dim} \left( x^i - \sum\limits_{j=1}^n x^i_j \phi_j(\xi,\eta,\chi) \right)^2 \\ \f$ 
     where \f$ x^i \f$ are the 3 coordinates (in 3D geometry case) of the node of interest.
     The function takes as input a double containing the node coordinates and a double 
     containing the \f$ (\xi,\eta,\chi) \f$ tuple.
     We adopt a linear approximation for the coordinates system, so the test functions 
     \f$ \phi_j(\xi,\eta,\chi) \f$ are the ones for linear interpolation.
     The function F is built in a way that the exact tuple \f$ (\tilde{\xi},\tilde{\eta},\tilde{\chi}) \f$, so that 
     \f$ \vec{x}(\tilde{\xi},\tilde{\eta},\tilde{\chi}) = \vec{x}_P \f$, represents an absolute minimum point of the function F.
     In particular \f$ F(\tilde{\xi},\tilde{\eta},\tilde{\chi}) = 0 \f$. 
     By setting a theresold value of F we decide when to stop the iterative calculation of the coordinates  \f$ (\xi,\eta,\chi) \f$.    
   */   
   double CalcF(double NodePos[],  /**< Coordinates of the target mesh node */
		double XiEtaChi[]  /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
               );
   
   //! Gradient of F 
   /*!
    * Function for the calculation of the gradient of F. The gradient components are written into the #_dF vector
    */
   void FirstDerF(double NodePos[],  /**< Coordinates of the target mesh node */
		  double XiEtaChi[]  /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
		 );
   
   //! Second derivatives of F 
   /*!
    * Function for the calculation of the Hessian matrix of F. The matrix components are written into the #_d2F vector
    */
   void SecondDerF(double NodePos[],  /**< Coordinates of the target mesh node */
		   double XiEtaChi[] /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
		  );
   
   
   //! Linear test function
   /*!
    * This function return the value of the 'node' linear test function on the point with coordinats 'XiEtaChi[]'
    */
//    double GenN(double XiEtaChi[], /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
// 	       int node      /**< Id of target mesh node */
// 	      );
   double LinPhi(int node,      /**< Id of target mesh node */
               double XiEtaChi[] /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */ 
              );
   //! First derivative of linear test function
   /*!
    * This function return the derivative value of the 'node' linear test function, on the point with coordinats 'XiEtaChi[]', in the direction 'direction'
    */
   double FirstDer(double XiEtaChi[], /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
		   int node,     /**< Id of target mesh node */
		   int direction /**< Direction along which the derivative is calculated */
		  );
   
   
   double SecondDer(double XiEtaChi[], /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
		    int node         /**< Id of target mesh node */
		   );                      
   
   //! Second derivative of linear test function
   /*!
    * This function Hessian matrix component H(row,col) of the linear test function 'node' on the node with coordinates 'XiEtaChi[]'
    */
   double SecondDer(double XiEtaChi[], /**< Array containing the \f$ (\xi,\eta,\chi) \f$ coordinates */
		    int node,          /**< Id of target mesh node */
		    int row,           /**< Row of Hessian matrix */
		    int col            /**< Column of the Hessian matrix */
		   );     
   
   
   //! Function that returns the number of source mesh nodes (per cell) used in the interpolation
   int  GetBoundNodesPerCell();
   
   //! Function that the domain type of the mesh group
   int  GetDomainType();
   
   //! Quadratic test function
   /*!
    * Function that returns the value of a particular test function a given point
    */
   double QuadPhi(int PhiNumber,  /**< Id of the test function */
		  double GPoint[] /**< Coordinates of the point in the reference element */
		 );
   
   
   //! Function for the interpolation
   /*!
    * This function calculates the interpolated field on the target mesh. 
    * The function return a MEDCoupling::MEDCouplingFieldDouble * object containing the interpolated solution for the target mesh
    */
   MEDCoupling::MEDCouplingFieldDouble * InterpolatedField(
      const MEDCoupling::MEDCouplingFieldDouble* SourceField, /**< Source mesh field containing the solution used for the interpolation */
      const int order = 2
   );
   
   // Funzione per interpolare con input aggiuntivo per il campo del problema FEMuS 
   
   MEDCoupling::MEDCouplingFieldDouble * InterpolatedField(
      const MEDCoupling::MEDCouplingFieldDouble* SourceField, /**< Source mesh field containing the solution used for the interpolation */
      const MEDCoupling::MEDCouplingFieldDouble* TargetField, /**< Target mesh field */
      const int order = 2
   );
   
   
    
   MEDCoupling::MEDCouplingFieldDouble * InterpolatedField(
      const MEDCoupling::MEDCouplingFieldDouble* SourceField, /**< Source mesh field containing the solution used for the interpolation */
      double DefaultValue,
      const int order = 2
   );
   
   
  /// \f$ \xi,\eta,\chi \f$ coordinates of the HEX8 vertices 
  const int _XiEtaChiVert[24] = {
    -1, 1, 1,-1,-1, 1, 1,-1,  /// xi
    -1,-1, 1, 1,-1,-1, 1, 1,  /// eta
    -1,-1,-1,-1, 1, 1, 1, 1   /// chi 
  }; 

  /// \f$ \xi,\eta,\chi \f$ coordinates of the HEX27 nodes
  const int _CoordHex27[27*3]={
  // 1   2    3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 27
    -1,  1,  1, -1, -1,  1,  1, -1,  0,  1,  0, -1,  0,  1,  0, -1, -1,  1,  1, -1,  0,  0,  1,  0, -1,  0, 0,  /// xi
    -1, -1,  1,  1, -1, -1,  1,  1, -1,  0,  1,  0, -1,  0,  1,  0, -1, -1,  1,  1,  0, -1,  0,  1,  0,  0, 0,  /// eta
    -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  1, 0   /// chi
  }; 
  
  /// \f$ \xi,\eta \f$ coordinates of the QUAD9 nodes
  const int _CoordQuad9[9*2]={
    /**< 1   2   3   4   5   6   7   8   9  node */
    -1,  1,  1, -1,  0,  1,  0, -1,  0, /**< coordinate xi */
    -1, -1,  1,  1, -1,  0,  1,  0,  0  /**< coordinate eta */
  }; 
  
  /// \f$ \xi \f$ coordinates of the EDGE3 nodes
  const int _CoordEdge3[3]={
  // 1   2   3   
    -1,  1,  0 /// xi
  }; 
   
  /// Offset for reading the generic node coordinates inside the _XiEtaChiVert array
  const int _VertCoordOff = 8; 
  
  /// Offset for reading the generic node coordinates inside the _CoordHex27 array
  const int _Hex27Off = 27; 
  
  /// Offset for reading the generic node coordinates inside the _CoordQuad9 array
  const int _Quad9Off = 9;  
  
  
  int _subsecAm;
  int _subsecAM;
  int _subsecBm;
  int _subsecBM;
  void SetSubSector(int sector, int iCell);
  short int* _Scount[4];
  short int* _SubSec[4];
  
 MEDCoupling::DataArrayDouble* _sec0;// = MEDCoupling::DataArrayDouble::New();
 MEDCoupling::DataArrayDouble* _sec2;// = MEDCoupling::DataArrayDouble::New();
 MEDCoupling::DataArrayDouble* _sec1;// = MEDCoupling::DataArrayDouble::New();
 MEDCoupling::DataArrayDouble* _sec3;// = MEDCoupling::DataArrayDouble::New();
  
 std::vector<double> _coord; 
 void CheckBelonging(bool &found);
 double _pos[3];
//   MEDCoupling::DataArrayDouble* _secArray;
 
  int IsFilled(){return __Filled;}
 
private:
  // Private parameters are characterised by double underscore
  const MEDCoupling::MEDCouplingUMesh * __InMesh;
  const MEDCoupling::MEDCouplingUMesh * __OutMesh;
  int  __BoundNodesPerCell;
  int  __TrgNodes;
  int  __Domain;
  int __Filled;
  
//   MEDCoupling::DataArrayInt * BoundingNodes;
//   MEDCoupling::DataArrayDouble * XiEta; 
};
#endif
 #endif
