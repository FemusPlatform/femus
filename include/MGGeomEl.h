#ifndef __mggeomel_h__
#define __mggeomel_h__

#include <string>
#include "MGFE_conf.h"
/// Class containing all the geometrical informations about the element.
/**Number of nodes and sides, both quadratic and linear, prolongation and embedding matrices
 and surface topology. */
class MGGeomEl  {

public:
  ///@{ \name GEOMETRIC CONFIGURATION OF THE ELEMENT 
     int n_q[2];             ///< Quadratic number of nodes (volume+surface)
     int n_l[2];             ///< Linear number of nodes (volume+surface)
     int n_se[2];            ///< Number of linear subelements (volume+surface)
     
     /// Number of element sides (volume+surface)
     int _n_sides[2];         
  ///@}
     
  ///@{ \name NAME OF THE ELEMENT
    std::string name[2];             ///< Element name (volume+surface)
    
    /// Print element name (volume+surface) for XDMF print (linear)
    std::string pname[2];            
  ///@} 
 
  ///@{ \name CONSTRUCTOR-DESTRUCTOR 
    MGGeomEl();  
    ~MGGeomEl();  
  ///@}

// ====== embedding matrices
#if ELTYPE==27  // Hex27 family ========================
#if DIMENSION==1
    static const double Prol[3*2];
    static const float _embedding_matrix_q[2][3][3];
    static const float _embedding_matrix_l[2][2][2];
    static const  int _surf_top[2];    ///< surface topology (side-> nodes)
#endif    
#if DIMENSION==2
    static const double Prol[9*4];
    static const float _embedding_matrix_q[4][9][9];
    static const float _embedding_matrix_l[4][4][4];
     static const  int _surf_top[12];    ///< surface topology (side-> nodes)
#endif     
#if DIMENSION==3    
    static const double Prol[27*8];
    static const float _embedding_matrix_q[8][27][27];
    static const float _embedding_matrix_l[8][8][8];// Hex8
    static const  int _surf_top[54];    ///< surface topology (side-> nodes)
#endif


#elif ELTYPE==10  // Tetra family ======================
#if DIMENSION==1
    static const double Prol[3*2];
    static const float _embedding_matrix_q[2][3][3];
    static const float _embedding_matrix_l[2][2][2];
    static const  int _surf_top[2];    ///< surface topology (side-> nodes)
#endif    
#if DIMENSION==2
    static const double Prol[6*3];
    static const float _embedding_matrix_q[4][6][6];// TRI6
    static const float _embedding_matrix_l[4][3][3];// TRI3
    static const  int _surf_top[9];    ///< surface topology (side-> nodes)
#endif    
#if DIMENSION==3
    static const double Prol[10*4];
    static const float _embedding_matrix_q[8][10][10];// TeT10
    static const float _embedding_matrix_l[8][4][4];// TET4
      static const  int _surf_top[24];    ///< surface topology (side-> nodes)
#endif


#elif ELTYPE==18  // Prism family ======================
#if DIMENSION==1
    static const double Prol[3*2];
    static const float _embedding_matrix_q[2][3][3];
    static const float _embedding_matrix_l[2][2][2];
    static const  int _surf_top[2];    ///< surface topology (side-> nodes)
#endif    
#if DIMENSION==2
    static const double Prol[6*3];
    static const float _embedding_matrix_q[4][6][6];// 
    static const float _embedding_matrix_l[4][3][3];// 
    static const  int _surf_top[9];    ///< surface topology (side-> nodes)
#endif    
#if DIMENSION==3
    static const double Prol[18*6];
    static const float _embedding_matrix_q[8][18][18];// 
    static const float _embedding_matrix_l[8][6][6];// 
    static const  int _surf_top[39];    ///< surface topology (side-> nodes)
#endif



#endif
};

#endif

