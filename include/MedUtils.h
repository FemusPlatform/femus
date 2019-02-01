#ifndef __MedUtilsg_h__
#define __MedUtilsg_h__


#include "Solverlib_conf.h"
#include <vector>
#include <map>

class MGFE;

#ifdef HAVE_MED
#include "MEDCouplingFieldDouble.hxx"

namespace MEDCoupling
{
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class DataArrayInt;
  class DataArrayDouble;
}

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

//! Method for Node to Cell conversion
/*!
 * This enumerator defines the criteria for the nodes to cell field conversion
 * in #GetCellField method.
 */
enum ConversionMode
{ MeanValue = 0, MidPoint = 1 };

/// Class for the interpolation of a solution from a source mesh to a target mesh
class MedUtils
{

public:

 int _proc;


  MedUtils ();

  //! Constructor of the MedUtils Class
  MedUtils (const MEDCoupling::MEDCouplingUMesh * SourceMesh, /**< Mesh support of the source geometry */
	const MEDCoupling::MEDCouplingUMesh * TargetMesh, /**< Mesh support of the target geometry */
	int DomainType = Boundary /**< Domain type of the mesh group (Boundary of Volume) */
    );
  //! Destructor of the MedUtils Class
//   ~MedUtils();
    virtual ~ MedUtils ();

  // PRINT FUNCTION
  void PrintMed (const MEDCoupling::MEDCouplingUMesh * SourceV_update,
		 std::vector < MEDCoupling::MEDCouplingFieldDouble * >f);

  void PrintMed (std::vector < MEDCoupling::MEDCouplingFieldDouble * >f, std::string FileName);

  void PrintMed (MEDCoupling::MEDCouplingFieldDouble * f,
		 std::string FileName, int n = 1);


  MGFE *_fe[3];
  void BuildCanonicalElementNodesMap (int NodesPerCell, std::map < int,
				      int >&Mappa);
  
 double Integrate (const MEDCoupling::MEDCouplingFieldDouble * Field, int order = 2, int n_cmp = 1,/**< Number of variables (in)   */
		    int first_cmp = 0,			     /**< First component     (in)   */
		    int method = 0,   /**< Method  #method (in) */
		    const MEDCoupling::MEDCouplingFieldDouble * VelField =  NULL);

  virtual void GaussLoop (
      std::vector < double >NodeVar,
			  std::vector < double >Velocity,
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
  MEDCoupling::MEDCouplingFieldDouble * GetCellField(const MEDCoupling::MEDCouplingFieldDouble* SourceField, int ConversionMode=MeanValue);
  void  InitFe();
  
  inline void setProcId(int procId){
      _proc = procId;
  }
  
// ======================================================================================================================  
  virtual MEDCoupling::MEDCouplingFieldDouble * InterpolatedField (
      const MEDCoupling::MEDCouplingFieldDouble * SourceField, /**< Source mesh field that have the solution used for the interpolation */
    int order = 2
);
// ======================================================================================================================
  virtual inline int IsFilled (){  return 0; };
 // ===================================================================================================================== 
  virtual void FillParameters (
      const MEDCoupling::MEDCouplingUMesh * SourceMesh, //< Mesh support of the source geometry 
      const MEDCoupling::MEDCouplingUMesh * TargetMesh, //< Mesh support of the target geometry 
      int DomainType = Boundary                         //< Domain type of the mesh group (Boundary of Volume) 
    );

};

#endif
#endif
