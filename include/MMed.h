#ifndef __MMedg_h__
#define __MMedg_h__


#include "Solverlib_conf.h"

#define USE_FEMUS (0)

#ifdef HAVE_MED

#include <vector>
#include<map>

// include local class


#ifdef USE_FEMUS==1
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGEquationsSystem.h"
#include "FEMUS.h"
#include "MGSolverBase.h"
#include "MeshExtended.h"
#include "MGTimeLoop.h"
#include "MGEquationsSystem.h"
#include "MGSystem.h"
#include "Equations_conf.h"
#include "MGGeomEl.h"
#include "MGFE_conf.h"
#endif


#include "MGFE.h"



// // // // // // // // // // // // // // // // // // // // // // // // // // // //

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}

// class FEMUS;
//! Method for the integration function #getMediumValuesOnBoundary_elem
enum method
{
  Mean = 0,    /**< Mean integral value */
  AxiMean = 1, /**< Mean integral value with axisymmetry */
  Bulk = 2,    /**< Mean integral value weighted on velocity component normal to surface */
  AxiBulk = 3, /**< Mean integral value weighted on velocity component normal to surface and axisymmetry */
  Area = 4,    /**< Area of the integration surface */
  Integrale = 5, /**< Integral value */
  AxiArea = 8,
  NormL2 = 16
};

//! Domain type of mesh group
/*!
 * Domain type of the mesh where the solution interpolation is performed.
 * We need to distinguish the interpolation on boundary groups and on volume groups because
 * the algorithm for searching the target node inside the source mesh is different
 */
enum DomainType
{ Boundary = 0, Volume = 1 };

/// Class for the interpolation of a solution from a source mesh to a target mesh
class MMed
{

public:




  MMed ();

  //! Constructor of the MMed Class
  MMed (const ParaMEDMEM::MEDCouplingUMesh * SourceMesh, /**< Mesh support of the source geometry */
	const ParaMEDMEM::MEDCouplingUMesh * TargetMesh, /**< Mesh support of the target geometry */
	int DomainType = Boundary /**< Domain type of the mesh group (Boundary of Volume) */
    );
  //! Destructor of the MMed Class
//   ~MMed();
    virtual ~ MMed ();

  // PRINT FUNCTION
  void PrintMed (const ParaMEDMEM::MEDCouplingUMesh * SourceV_update,
		 std::vector < ParaMEDMEM::MEDCouplingFieldDouble * >f);
#ifdef USE_FEMUS==1
  void PrintMed (FEMUS * PFemus, int n = 0);
  void PrintMed (std::vector < FEMUS * >PFemus);
  void PrintMed (std::vector < ParaMEDMEM::MEDCouplingFieldDouble * >f,
		 std::string FileName);
#endif
  void PrintMed (ParaMEDMEM::MEDCouplingFieldDouble * f,
		 std::string FileName, int n = 1);


  MGFE *_fe[3];
  void BuildCanonicalElementNodesMap (int NodesPerCell, std::map < int,
				      int >&Mappa);

#ifdef USE_FEMUS==1
  virtual double Integrate (FEMUS * PFemus,  /**< FEMus problem           (in)*/
			    int id,	   /**< Interface name      (in)   */
			    const char *system_name,
					   /**< Equation name       (in)   */
			    int n_cmp,	   /**< Number of variables (in)   */
			    int first_cmp = 0,
					   /**< First component     (in)   */
			    int method = 0 /**< Method  #method (in) */
    );
#endif
  double Integrate (const ParaMEDMEM::MEDCouplingFieldDouble * Field, int order = 2, int n_cmp = 1,
				       /**< Number of variables (in)   */
		    int first_cmp = 0,
				     /**< First component     (in)   */
		    int method = 0,   /**< Method  #method (in) */
		    const ParaMEDMEM::MEDCouplingFieldDouble * VelField =
		    NULL);

#ifdef USE_FEMUS==1
  virtual ParaMEDMEM::MEDCouplingFieldDouble * GetVelocityField (FEMUS * PFem,	///[in] Pointer to FEMuS problem
								 int InterfaceId,	///[in] Interface ID from where we get the velocity field
								 int IsVelCoupled	///[in] Flag (1 or 0) to say if the velocity field is calculated with coupled (1) or uncoupled (0) system
    );
#endif
  virtual void GaussLoop (vector < double >NodeVar,
			  vector < double >Velocity,
			  bool BulkMedium,
			  const int dim,
			  const int DimRelToMax,
			  const int XOrd,
			  const int FOrd,
			  const int VOrd,
			  const int NodesPerCell,
			  std::vector < double >PointsCoords, int Fcc);


  int _CylCoord;
  int _BulkMedium;
  int _AreaCalc;
  int _IntegCalc;
  int _L2Norm;
  double _INTEGRAL, _AREA, _VELOCITY;
  void IntCoefficients (int rad);
  void FieldNodes (const int NodesPerCell, int &Fcc, int order);
  double GetIntResult (int method);

  void InitFe ();
  virtual ParaMEDMEM::MEDCouplingFieldDouble * InterpolatedField (const ParaMEDMEM::MEDCouplingFieldDouble * SourceField,
							       /**< Source mesh field containing the solution used for the interpolation */
								  int order =
								  2);

  virtual inline bool IsFilled ()
  {
    return false;
  };
  virtual void FillParameters (const ParaMEDMEM::MEDCouplingUMesh * SourceMesh,
							 /**< Mesh support of the source geometry */
			       const ParaMEDMEM::MEDCouplingUMesh * TargetMesh,
							 /**< Mesh support of the target geometry */
			       int DomainType = Boundary
				  /**< Domain type of the mesh group (Boundary of Volume) */
    );
  #ifdef USE_FEMUS==1
  void CreateInterfaces (FEMUS * PFemus);
#endif
};
#endif
#endif
