#ifndef __IbUtils__
#define __IbUtils__

#include "Domain_conf.h"
#include "InterfaceProjection.h"
#include "InterpolationOptions.hxx"

#include <map>
#include <sstream>
#include <vector>

namespace MEDCoupling {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
class MEDCouplingRemapper;
}  // namespace MEDCoupling

class IbUtils : public InterfaceProjection {
 private:
  const int _Proc;
  int _LibToMed_3D[27] = {7, 4,  5,  6,  3,  0,  1,  2,  19, 16, 17, 18, 11, 8,
                          9, 10, 15, 12, 13, 14, 25, 24, 21, 22, 23, 20, 26};
  int _LibToMed_2D[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  int** _CanElemMap;
  //   int CanonicalElementMaps
  MEDCoupling::MEDCouplingRemapper* _PieceWiseInterpolator;

 public:
  INTERP_KERNEL::IntersectionType _PieceWiseInterpScheme;
  MEDCoupling::MEDCouplingUMesh* _SolidBodyMesh;
  const MEDCoupling::MEDCouplingUMesh* _FluidMesh;

  // linear meshes for 3d med routines
  MEDCoupling::MEDCouplingUMesh* _LinearSolidBodyMesh;
  MEDCoupling::MEDCouplingUMesh* _LinearFluidMesh;

  MEDCoupling::MEDCouplingUMesh* _ExtractedMesh;
  MEDCoupling::MEDCouplingUMesh* _ReconInterface;
  MEDCoupling::MEDCouplingUMesh* _TriangulatedNearBoundary;

  MEDCoupling::MEDCouplingFieldDouble* _SolidCellColor;
  MEDCoupling::MEDCouplingFieldDouble* _SolidNodeColor;
  MEDCoupling::MEDCouplingFieldDouble* _SolidVelocity;

  MEDCoupling::MEDCouplingFieldDouble* _FluidCellColor;
  MEDCoupling::MEDCouplingFieldDouble* _FluidNodeColor;
  MEDCoupling::MEDCouplingFieldDouble* _ExtractedColor;
  MEDCoupling::MEDCouplingFieldDouble* _ExtractedVelocity;
  MEDCoupling::MEDCouplingFieldDouble* _InterpolatedSolidVelocity;
  MEDCoupling::MEDCouplingFieldDouble* _TriangulatedNodesIds;

  int* _CellToExtract;
  int _NumCellToExtract;
  double* _InterfaceNormal;
  double _Velocity[DIMENSION];
  double _SolidArea;
  int _iter;
  std::map<std::string, std::string> _FileMap;  /// String map containing Tproperties.in parameters
  double _rhof, _rhos, _gravity[3], _muf;
  int _buoyancy;
  std::map<int, std::vector<int>> _TriangulatedCellsPerSourceCell;
  std::map<int, std::vector<int>> _SolidIds;

  double _LowerThreshold;
  double _UpperThreshold;
  int _SimpleIBMethod;
  std::vector<std::map<int, double>> _RemapMatrix;

 public:
  IbUtils();
  IbUtils(
      MEDCoupling::MEDCouplingUMesh* SolidBodyMesh, const MEDCoupling::MEDCouplingUMesh* FluidMesh,
      int proc = 0);
  IbUtils(
      MEDCoupling::MEDCouplingUMesh* SolidBodyMesh, const MEDCoupling::MEDCouplingUMesh* FluidMesh,
      int p_femus, int proc);
  ~IbUtils();
  void read_file(int p_femus = 0);
  void print_par();
  void SetMeshes(
      MEDCoupling::MEDCouplingUMesh* SolidBodyMesh, const MEDCoupling::MEDCouplingUMesh* FluidMesh);

  void MoveMeshWithStress(MEDCoupling::MEDCouplingFieldDouble* VelocityField, double dt);
  void ComputeStress(MEDCoupling::MEDCouplingFieldDouble* ReducedVelocityField, double Stress[]);

  MEDCoupling::MEDCouplingFieldDouble* ComputeStressField(
      MEDCoupling::MEDCouplingFieldDouble* ReducedVelocityField,
      MEDCoupling::MEDCouplingFieldDouble* BoundaryField, const MEDCoupling::MEDCouplingUMesh* BoundaryMesh,
      const int axisym = 0);

  void InitColor();
  void CalcInterface(MEDCoupling::MEDCouplingFieldDouble* ProjectedColor);
  void CalcInterface2(MEDCoupling::MEDCouplingFieldDouble* ProjectedColor);

  void InitSolidVelocity(std::vector<double> Velocity);
  void InitPosition(double FirstTranslation[], MEDCoupling::MEDCouplingFieldDouble* VelocityField);
  void MoveWithGivenVectorAndRotation(
      double AngVel, double center[], double vector[], double dt, double velocity[],
      MEDCoupling::MEDCouplingFieldDouble* VelocityField);

  void UpdateFluidMesh(const MEDCoupling::MEDCouplingUMesh* NewMesh);
  void UpdateSolidMesh(const MEDCoupling::MEDCouplingUMesh* NewMesh);
  void UpdateSolidVelocity(const MEDCoupling::MEDCouplingFieldDouble* SourceField);
  void UpdateInterpolatedFields(MEDCoupling::MEDCouplingFieldDouble* VelocityField);
  void InsertTriangulatedCell(
      int SourceCellId, std::vector<int>& NodesPerCell, std::vector<int> OrderedNodes,
      std::vector<int>& Sequence, std::vector<double>& TriangulationCoords, double FullCoordinates[],
      int offset);

  MEDCoupling::MEDCouplingFieldDouble* InterpolateFluidOnSolid(
      MEDCoupling::MEDCouplingFieldDouble* VelocityField,
      MEDCoupling::MEDCouplingFieldDouble* SolidBoundaryField,
      const MEDCoupling::MEDCouplingUMesh* BoundaryMesh);
  MEDCoupling::MEDCouplingFieldDouble* InterpolateSolidOnFluid(
      MEDCoupling::MEDCouplingFieldDouble* SolidVelocity, MEDCoupling::MEDCouplingFieldDouble* FluidVelocity);
  MEDCoupling::MEDCouplingFieldDouble* GetFluidCellColor();
  MEDCoupling::MEDCouplingFieldDouble* GetFluidNodeColor();
  MEDCoupling::MEDCouplingFieldDouble* GetInterpolatedSolidVelocity();

  void GetCellContainingPoint(int TargetCellId, double NodeCoord[], double CanPos[], int& CellId);

  double TriangleArea(double* point1, double* point2, double* point3);
  double TriangleArea(double* point1);
  double QuadrangleArea(double* point1);
};

#endif
