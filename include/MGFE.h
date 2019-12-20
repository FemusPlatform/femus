#ifndef __mggauss_h__
#define __mggauss_h__

#include        "Domain_conf.h"

#include "Printinfo_conf.h"

#include <fstream>
#include <vector>
#include "MGUtils.h" 


/// Class containing mathematical information about the finite element.
/** 
 * 
 * */
class MGFE
#ifdef LM_REFCOUNT
      : public ReferenceCountedObject<MGFE>
#endif
{

public:

//   const MGUtils& _mgutils; ///< MGUtils pointer

  ///@{ \name FINITE ELEMENT DEFINITION 
  
  const int _dim; ///< Dimension (3 or 2) used in the 3D FEM
  
  /// Order of the shape functions (2=phi=QUADRATIC, 1=psi=LINEAR)
  const int _order;
  ///@}
  int _FamilyType;

  ///@{ \name FEM INTEGRATION 
  int _NoGauss1[3]; ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  /// Number of shape functions in 1-2-3D
  int _NoShape[3]; 
  /// Polinomial degree
  int _deg[3];
  
  // MODIFICA -> AGGIUNTA PER ELEMENTI MISTI TRIANGOLI - RETTANGOLI
  int _GNoGauss1[2][3]; ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  /// Number of shape functions in 1-2-3D
  int _GNoShape[2][3]; 
  /// Polinomial degree
  int _Gdeg[2][3];  
  
  
  ///@}
  
  ///@{ \name  SHAPES [0]=1D,[1]=2D,[2]=3D in Gaussian points
  double* _weight1[3];              ///< Weight
  double* _phi_map1[3];             ///< Shape functions
  double* _dphidxez_map1[3];        ///< Shape derivative functions in gaussian points
  double* _dphidxez_map1_nodes[3];  ///< Shape derivative functions in nodal points
  double* _dphidxx_map1[3];         ///< Second order shape derivatives in Gaussian points
  
  
  double* _Gweight1[6];              ///< Weight
  double* _Gphi_map1[6];             ///< Shape functions
  double* _Gdphidxez_map1[6];        ///< Shape derivative functions in gaussian points
  double* _Gdphidxez_map1_nodes[6];  ///< Shape derivative functions in nodal points
  double* _Gdphidxx_map1[6];         ///< Second order shape derivatives in Gaussian points  
  int _GlobalFE;
  ///@}

//------------------------------------------------------------
  ///@{ \name CONSTRUCTOR-DESTRUCTOR 
  MGFE(
//    const MGUtils& mgutils_in,  ///< MGUtils pointer
   const int order,          ///< Order of the shape functions
   
   /// Fem type: 27(Quad,Hex) or 10(Tri,Tet)
   const int fem_type         
  );  
  ~MGFE();
  ///@}
  
  MGFE(const int order); 
  
  
  ///@{ \name INITIALIZING AND CLEANING
  void init_qua();    ///< Generates the Lagrangian quad shape functions
  void init_qua_tri(); 
  void init_qua_Gtri(); 
  void init_qua_rec(); 
  void init_qua_Grec();
  
  void init_lin();    ///< Generates the Lagrangian linear shape functions  
  void init_lin_tri();
  void init_lin_Gtri();
  void init_lin_rec();
  void init_lin_Grec();
  
  // -------------------------------------------------------
  // Triangle - based elements: 
  
  double Tri_2d_LinearPhi(int nPhi, double point[]);
  double Tri_2d_LinearDerPhi(int nPhi, double point[], int dir);
  double Tri_2d_QuadraticPhi(int nPhi, double point[]);
  double Tri_2d_QuadraticDerPhi(int nPhi, double point[], int dir);
  double Tri_2d_QuadraticDer2Phi(int nPhi, double point[], int dir1, int dir2);
  
  double Tri_3d_LinearPhi ( int nPhi, double point[] ); 
  double Tri_3d_LinearDerPhi ( int nPhi, double point[], int dir ) ;
  double Tri_3d_QuadraticPhi ( int nPhi, double point[] );
  double Tri_3d_QuadraticDerPhi ( int nPhi, double point[], int dir );
  double Tri_3d_QuadraticDer2Phi ( int nPhi, double point[], int dir1, int dir2 );
  
  // -------------------------------------------------------
  // Quadrangle - based elements: 
  double Edge_Lin_Phi(int PhiCoeff, double Coordinate);
  double Edge_Lin_DPhi(int PhiCoeff, double Coordinate);
  double Edge_Quad_Phi(int PhiCoeff, double Coordinate);
  double Edge_Quad_DPhi(int PhiCoeff, double Coordinate);
  double Edge_Quad_D2Phi(int PhiCoeff, double Coordinate);
  
  double Rec_Lin_Phi(int nPhi, double point[], int dimension);
  double Rec_Lin_DPhi(int nPhi, double point[], int dimension, int DirDer);
  double Rec_Lin_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2);
  double Rec_Quad_Phi(int nPhi, double point[], int dimension);
  double Rec_Quad_DPhi(int nPhi, double point[], int dimension, int DirDer);
  double Rec_Quad_D2Phi(int nPhi, double point[], int dimension, int DirDer1, int DirDer2);
  
  double FirstDerivateOfLocalPhi(int nPhi, double point[], int dimension, int DirDer, int FamilyType);
  int GetFamilyType(int elem_dof, int dim); // the function returns 0 or 1
  
  void init_pie();    ///< Generates the Lagrangian piecewise shape functions
  
  const int _MedToLib_27[27]= {
        //0  1   2   3   4   5   6   7   8    9    10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 
          0, 3,  2,  1,  4,  7,  6,  5,  11,  10,  9,  8,  19, 18, 17, 16, 12, 15, 14, 13, 20, 24, 23, 22, 21, 25, 26
    };
    
  
    /// \f$ \xi,\eta,\chi \f$ coordinates of the HEX27 nodes
    const int _CooH27[27*3]= {
        //0  1   2    3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 
        -1,  1,  1, -1, -1,  1,  1, -1,  0,  1,  0, -1, -1,  1,  1, -1,  0,  1,  0, -1,  0,  0,  1,  0, -1,  0, 0,  /// xi
        -1, -1,  1,  1, -1, -1,  1,  1, -1,  0,  1,  0, -1, -1,  1,  1, -1,  0,  1,  0,  0, -1,  0,  1,  0,  0, 0,  /// eta
        -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,  1,  1,  1,  1, -1,  0,  0,  0,  0,  1, 0   /// chi
    };
    /// \f$ \xi,\eta \f$ coordinates of the QUAD9 nodes
    const int _CooQ9[9*2]= {
        /**< 1   2   3   4   5   6   7   8   9  node */
        -1,  1,  1, -1,  0,  1,  0, -1,  0, /**< coordinate xi */
        -1, -1,  1,  1, -1,  0,  1,  0,  0  /**< coordinate eta */
    };
    /// \f$ \xi \f$ coordinates of the EDGE3 nodes
    const int _CooE3[3]= {
      // 1   2   3
        -1,  1,  0 /// xi
    };
    
    const int _CooTriEl[9] = {// L_i = a_i*xi + b_i*eta + c_i
  //   a   b   c
      -1, -1,  1,  // 0
       1,  0,  0,  // 1
       0,  1,  0   // 2
    };
    
    const int _CooTetra4[16] = {// L_i = a_i*xi + b_i*eta + c_i*chi + d_i
  //   a   b   c   d
      -1, -1, -1,  1, // 0
       1,  0,  0,  0, // 1
       0,  1,  0,  0, // 2
       0,  0,  1,  0  // 3
    };
    
    const int _CooTetra10[40] = {// L_i = l(j)*(a* l(k) + b)
  //   j   k   a    b
       0,  0,  2,  -1, // 0
       1,  1,  2,  -1, // 1
       2,  2,  2,  -1, // 2
       3,  3,  2,  -1, // 3
       1,  0,  4,   0, // 4
       1,  2,  4,   0, // 5
       2,  0,  4,   0, // 6
       0,  3,  4,   0, // 7
       1,  3,  4,   0, // 8
       2,  3,  4,   0  // 9
    }; 
    /// Offset for reading the generic node coordinates inside the _CooH27 array
    const int _H27Off = 27;
    /// Offset for reading the generic node coordinates inside the _CooQ9 array
    const int _Q9Off = 9;
  /// Clear data substructures
  void clear();  
  ///@}

  
  ///@{ \NAME RETURN SHAPE DERIVATIVE FUNCTIONS IN GAUSS POINTS
   ///< \param[in]  <dim>    Dimension
   ///< \param[in]  <qp>     Gaussian point
   ///< \param[in]  <InvJac> Jacobean
   ///< \param[out] <dphi>   Derivative
  
  void get_dphi_gl_g(const int dim,const int qp,const  double InvJac[], double dphi[]);
  void get_ddphi_gl_g(const int dim,const int qp,const double InvJac[], double ddphi[]);
  void get_dphi_gl_g(const int dim,const int qp,const  double InvJac[], double dphi[],int sdim);
  void get_dphi_gl_g(const int dim,const int qp, const double InvJac[],std::vector<double> &dphi);
  void get_dphi_node(const int dim,const int qp,const  double InvJac[], double dphi[]);
  
  void get_dphi_on_given_node(const int dim, double ElemCoords[], double CanPos[], double dphi[]);
  void get_dphi_on_given_nodeG(const int dim, double ElemCoords[], double CanPos[], double dphi[], int FamilyType);
  
  void get_dphi_arb_node(std::vector<double> NodeCoord, const int order, double InvJac[], double dphi[]);
  ///@}
  
  ///@{ \NAME RETURN SHAPE FUNCTIONS IN GAUSS POINTS
   ///< \param[in]  <dim>    Dimension
   ///< \param[in]  <qp>     Gaussian point
   ///< \param[out] <phi>    Shape function
  
  void get_phi_gl_g(const int dim, const int qp,double phi[]);
  void get_phi_gl_g(const int dim, const int qp, std::vector<double>& phi);
  ///@}
      
  void get_phi_g_arb_el(const int dim, const int qp,double phi[], int FamilyType);
  void get_dphi_g_arb_el(const int dim,const int qp,const  double InvJac[], double dphi[], int FamilyType);
  
  //  Jacobian functions ---------------------------------
  double ComputeInverseMatrix(double Matrix[], double InvMatrix[], int Dimension);
  // Jacobian at gaussian points
  double Jac(const int ng, double x[], double InvJac[]) ;  ///< Jacobian (3D-2D)
  double Jac1D(const int ng,double x[], double InvJac[]);///< Jacobian (1D)
  double Jac2D(const int ng,double x[], double InvJac[]);  ///< Jacobian (2D)
  double Jac2D(const int ng,double x[], double InvJac[], double Jac[]);  ///< Jacobian (2D)
  double Jac3D(const int ng,double x[], double InvJac[]);///< Jacobian (3D)
  
  // Jacobian at gaussian points for arbitrary fem element
  double JacG(const int ng, double x[], double InvJac[], int FamilyType, int dimension) ;  ///< Jacobian (3D-2D)
  
  // Jacobian at nodal points
  double Jac_nodes(const int ng, double x[], double InvJac[]) ;  ///< Jacobian (3D-2D)
  double Jac1D_nodes(const int ng,double x[], double InvJac[]);///< Jacobian (1D)
  double Jac2D_nodes(const int ng,double x[], double InvJac[]);  ///< Jacobian (2D)
  double Jac3D_nodes(const int ng,double x[], double InvJac[]);///< Jacobian (3D)

  void JacOnGivenCanCoords(int dim, double ElemCoords[], double CanCoords[], double InvJac[], int FamilyType, int nShape);

  // Jacobian dim -1
  double JacSur(const int ng, double x[], double InvJac[]) const;      ///< Boundary Jacobian (2D-1D)
  double JacSur2D(const int ng, double x[], double InvJac[]) const;
  double JacSur3D(const int n_gauss,double x[], double InvJac[]) const;///< Surface Jacobian (2D)
//   double JacSur2D(const int ng,double x[]) const;     ///<  Line Jacobian    (1D)
  double JacSur1D(const int ng,double x[], double InvJac[]) const;     ///<  P Jacobian    (1D)
  // functions -----------------------------------------
  /// Compute normal normal_g[] at xx[] point
  void normal_g(const double* xx, double* normal_g) const; ///< unit normal to the surface
 void normal_g(const double* xx,const double x_c[],double* normal_g) const;///<  and check x_c[] being interior point
  void normal_g(const double* xx,const double x_c[],double* normal_g, int & sign_normal) const;///<  and check x_c[] being interior point
  // Reading - writing -------------------------------------
  /// Write
  void write(const std::string& name,const int kdim);
private:
  /// Reading
  void read_c(std::istream& infile, const int kdim);
  /// Writing
  void write_c(std::ostream& infile, const int kdim);
};

#endif
