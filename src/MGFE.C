// std library -----------------
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// class headers ----------
#include "MGFE.h"           // class header
#include "MGFE_conf.h"       // configuration

// local ------------------
#include "MGUtils.h"        // data and path to file
#include "Printinfo_conf.h" // to print info


// ==============================
/// This function constructs the class MGFE:
MGFE::MGFE (
    const int order_in,        // Fem order input
    const int fem_type         // Fem type: 27(Quad,Hex) or 10(Tri,Tet)
) :
_dim ( DIMENSION ),      // space dimension
_order ( order_in ) {    // Fem order
    // =============================
    int n_shape[3];
    int deg[3];
    int n_gauss[3];
    _GlobalFE = 0;
    _FamilyType = ( fem_type==27 ) ? 1:0;
    if ( order_in == 2 )
        switch ( fem_type ) {
        case 27:
            n_shape[0]=3;
            n_shape[1]=9;
            n_shape[2]=27;   // hex quadratic  shapes
            deg[0]=2;
            deg[1]=4;
            deg[2]=8;        // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=9;
            n_gauss[2]=27;   // 3x3x3 gaussian points
            break;
        case 10:
            n_shape[0]=3;
            n_shape[1]=6;
            n_shape[2]=10;   // tet quadratic  shapes
            deg[0]=2;
            deg[1]=2;
            deg[2]=3;        // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=7;
            n_gauss[2]=14;   // 3x3x3 gaussian points
            break;
        default :
            std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
            break;
        }
    else if ( order_in == 1 )
        switch ( fem_type ) {
        case 27:
            n_shape[0]=2;
            n_shape[1]=4;
            n_shape[2]=8;   // hex linear  shapes
            deg[0]=1;
            deg[1]=2;
            deg[2]=3;        // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=9;
            n_gauss[2]=27;   // 3x3x3 gaussian points
            break;
        case 10:
            n_shape[0]=2;
            n_shape[1]=3;
            n_shape[2]=4;   // tet linear  shapes
            deg[0]=1;
            deg[1]=2;
            deg[2]=3;        // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=7;
            n_gauss[2]=14;   // 3x3x3 gaussian points
            break;
        default :
            std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
            break;
        }
    else if ( order_in == 0 )
        switch ( fem_type ) {
        case 27:
            n_shape[0]= ( ( NDOF_K>1 ) ? 2  : 1 );
            n_shape[1]= ( ( NDOF_K>1 ) ?3 : 1 );
            n_shape[2]= ( ( NDOF_K>1 ) ? 4  : 1 ); // hex piecewise shapes
            deg[0]= ( ( NDOF_K>0 ) ? 1 : 0 );
            deg[1]= ( ( NDOF_K>0 ) ? 1 : 0 );
            deg[2]= ( ( NDOF_K>0 ) ? 1 : 0 );  // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=9;
            n_gauss[2]=27;   // 3x3x3 gaussian points
            break;
        case 10:
            n_shape[0]=1;
            n_shape[1]=1;
            n_shape[2]=1;   // tet piecewise shapes
            deg[0]=0;
            deg[1]=0;
            deg[2]=0;        // degree of the shape
            n_gauss[0]=3;
            n_gauss[1]=7;
            n_gauss[2]=14;   // 3x3x3 gaussian points
            break;
        default :
            std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";
            break;
        }
    else std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";

// Fem -------------------------
    for ( int idim=0; idim<3; idim++ ) {
        _NoGauss1[idim]=n_gauss[idim];// # of Gaussian points
        _NoShape[idim]=n_shape[idim]; //  # of shape functions
        _deg[idim]=deg[idim];
    }

// weights, shapes and shape derivatives -----------------------------
    for ( int idim=0; idim<3; idim++ ) {
        _weight1[idim]= new double [_NoGauss1[idim]];
        _phi_map1[idim]=new double[_NoShape[idim]*_NoGauss1[idim]];
        _dphidxez_map1[idim]=new double[ ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
        _dphidxez_map1_nodes[idim]=new double[ ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
        _dphidxx_map1[idim]=new double[ ( idim+1 ) * ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
    }

#ifdef PRINT_INFO
    std::cout << "\n MGFE::MGFE: gaussian points: "
    << _NoGauss1[DIMENSION-1] << " \n";
    std::cout << " MGFE::MGFE: polynomial order " << _order
    << " - # shape func=(";
    for ( int idim=0; idim<3; idim++ ) std::cout <<	_NoShape[idim] << ",";
    std::cout    << ") \n";
#endif
    return;
}

MGFE::MGFE (
    const int order_in        // Fem order input
) :
_dim ( DIMENSION ),      // space dimension
_order ( order_in ) {    // Fem order
    // =============================

    int n_shape[3];
    int deg[3];
    int n_gauss[3];
    _GlobalFE = 1;

// _GNoShape[x][i], _Gdeg[x][i], _GNoGauss1[x][i]  ->  x is 0 (tri-based) or 1 (rect-based)
    if ( order_in == 2 ) {
// HEXA
        _GNoShape[1][0]=3;
        _GNoShape[1][1]=9;
        _GNoShape[1][2]=27;   // hex quadratic  shapes
        _Gdeg[1][0]=2;
        _Gdeg[1][1]=4;
        _Gdeg[1][2]=8;        // degree of the shape
        _GNoGauss1[1][0]=3;
        _GNoGauss1[1][1]=9;
        _GNoGauss1[1][2]=27;   // 3x3x3 gaussian points
// TRI
        _GNoShape[0][0]=3;
        _GNoShape[0][1]=6;
        _GNoShape[0][2]=10;   // tet quadratic  shapes
        _Gdeg[0][0]=2;
        _Gdeg[0][1]=2;
        _Gdeg[0][2]=3;        // degree of the shape
        _GNoGauss1[0][0]=3;
        _GNoGauss1[0][1]=7;
        _GNoGauss1[0][2]=14;   // 3x3x3 gaussian points

    } else if ( order_in == 1 ) {
// HEXA
        _GNoShape[1][0]=2;
        _GNoShape[1][1]=4;
        _GNoShape[1][2]=8;   // hex linear  shapes
        _Gdeg[1][0]=1;
        _Gdeg[1][1]=2;
        _Gdeg[1][2]=3;        // degree of the shape
        _GNoGauss1[1][0]=3;
        _GNoGauss1[1][1]=9;
        _GNoGauss1[1][2]=27;   // 3x3x3 gaussian points
//TRI
        _GNoShape[0][0]=2;
        _GNoShape[0][1]=3;
        _GNoShape[0][2]=4;   // tet linear  shapes
        _Gdeg[0][0]=1;
        _Gdeg[0][1]=2;
        _Gdeg[0][2]=3;        // degree of the shape
        _GNoGauss1[0][0]=3;
        _GNoGauss1[0][1]=7;
        _GNoGauss1[0][2]=14;   // 3x3x3 gaussian points
    } else if ( order_in == 0 ) {
// HEXA
        _GNoShape[1][0]= ( ( NDOF_K>1 ) ? 2  : 1 );
        _GNoShape[1][1]= ( ( NDOF_K>1 ) ?3 : 1 );
        _GNoShape[1][2]= ( ( NDOF_K>1 ) ? 4  : 1 ); // hex piecewise shapes
        _Gdeg[1][0]= ( ( NDOF_K>0 ) ? 1 : 0 );
        _Gdeg[1][1]= ( ( NDOF_K>0 ) ? 1 : 0 );
        _Gdeg[1][2]= ( ( NDOF_K>0 ) ? 1 : 0 );  // degree of the shape
        _GNoGauss1[1][0]=3;
        _GNoGauss1[1][1]=9;
        _GNoGauss1[1][2]=27;   // 3x3x3 gaussian points
//TRI
        _GNoShape[0][0]=1;
        _GNoShape[0][1]=1;
        _GNoShape[0][2]=1;   // tet piecewise shapes
        _Gdeg[0][0]=0;
        _Gdeg[0][1]=0;
        _Gdeg[0][2]=0;        // degree of the shape
        _GNoGauss1[0][0]=3;
        _GNoGauss1[0][1]=7;
        _GNoGauss1[0][2]=14;   // 3x3x3 gaussian points

    } else std::cout << "ERROR! FEM SHAPES BAD DEFINED"<< "\n";

// Fem -------------------------
    for ( int idim=0; idim<3; idim++ ) {
        _NoGauss1[idim]= _GNoGauss1[1][idim];
        _NoShape[idim]=_GNoShape[1][idim];
        _deg[idim]=_Gdeg[1][idim];
    }
// weights, shapes and shape derivatives -----------------------------
    for ( int idim=0; idim<3; idim++ ) {
        _weight1[idim]             = new double [_NoGauss1[idim]];
        _phi_map1[idim]            = new double[_NoShape[idim]*_NoGauss1[idim]];
        _dphidxez_map1[idim]       = new double[ ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
        _dphidxez_map1_nodes[idim] = new double[ ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
        _dphidxx_map1[idim]        = new double[ ( idim+1 ) * ( idim+1 ) *_NoShape[idim]*_NoGauss1[idim]];
    }
    // global structures - triangle and rectangle based finite element
    for ( int i_type = 0; i_type<2; i_type++ ) {
        for ( int idim=0; idim<3; idim++ ) {
            int comp = idim + i_type*3;
            _Gweight1[comp]             = new double[_GNoGauss1[i_type][idim]];
            _Gphi_map1[comp]            = new double[_GNoShape[i_type][idim]*_GNoGauss1[i_type][idim]];
            _Gdphidxez_map1[comp]       = new double[ ( idim+1 ) *_GNoShape[i_type][idim]*_GNoGauss1[i_type][idim]];
            _Gdphidxez_map1_nodes[comp] = new double[ ( idim+1 ) *_GNoShape[i_type][idim]*_GNoGauss1[i_type][idim]];
            _Gdphidxx_map1[comp]        = new double[ ( idim+1 ) * ( idim+1 ) *_GNoShape[i_type][idim]*_GNoGauss1[i_type][idim]];
        }
    }
#ifdef PRINT_INFO
    std::cout << "\n MGFE::MGFE: gaussian points: triangle-based " << _GNoGauss1[0][DIMENSION-1] <<" rectangle-based "<<_GNoGauss1[1][DIMENSION-1]<< " \n";
    std::cout << " MGFE::MGFE: polynomial order " << _order << " - # triangle-based FE shape func=(";
    for ( int idim=0; idim<3; idim++ ) std::cout <<	_GNoShape[0][idim] << ",";
    std::cout    << ") \n";
    std::cout << " MGFE::MGFE: polynomial order " << _order << " - # rectangle-based FE shape func=(";
    for ( int idim=0; idim<3; idim++ ) std::cout <<	_GNoShape[1][idim] << ",";
    std::cout    << ") \n";
#endif
    return;
}
// ==============================
/// This function destroys the MGFE class
MGFE::~MGFE (
) {// ==============================
    clear();
}

// ================================
/// This function clears the substructures
void MGFE::clear (
) {// ========================
    for ( int idim=0; idim<_dim; idim++ ) {
        delete []_weight1[idim];       // weight
        delete []_phi_map1[idim];      // shapes at g.p.
        delete []_dphidxez_map1[idim]; // derivatives at g.p.
        delete []_dphidxez_map1_nodes[idim]; // derivatives on nodes
        delete []_dphidxx_map1[idim]; // derivatives on nodes
        if ( _GlobalFE==1 ) {
            delete [] _Gdphidxx_map1[idim];
            delete [] _Gdphidxx_map1[idim + 3];
            delete [] _Gweight1[idim];
            delete [] _Gweight1[idim + 3];
            delete [] _Gphi_map1[idim];
            delete [] _Gdphidxez_map1[idim];
            delete [] _Gdphidxez_map1_nodes[idim];
            delete [] _Gphi_map1[idim +3];
            delete [] _Gdphidxez_map1[idim +3];
            delete [] _Gdphidxez_map1_nodes[idim +3];
        }
    }
    return;
}


/// /// This function generates the Lagrangian quad shape functions
void MGFE::init_qua (
) {// ================================
    if ( _GlobalFE==1 ) {
        init_qua_Grec();
        init_qua_Gtri();
    }
//     else {
#if ELTYPE==27
    init_qua_rec();
#endif

#if ELTYPE==10
    init_qua_tri();
#endif
//     }
    return;
}

void MGFE::init_qua_rec() {
//               ********************************************
//                               QUAD 9
// 		 ********************************************
//
// 			       3 ______6_____ 2
// 				|            |
// 			        |            |
//			       7|      8     |5
//				|            |
// 			        |____________|
//			       0       4      1

//gaussian coordinates
    const double a=-sqrt ( 3./5. );
    const double b=0.;
    const double c=-a;
    const double a_n=-1.;
    const double b_n=0.;
    const double c_n=1.;

    const double x1D[3]= {a,b,c};
    const double x1D_n[3]= {a_n,b_n,c_n};

    const double x2D[9]= {a,a,a,b,b,b,c,c,c};
    const double y2D[9]= {a,b,c,a,b,c,a,b,c};
    const double x2D_n[9]= {a_n,c_n,c_n,a_n,b_n,c_n,b_n,a_n,b_n};
    const double y2D_n[9]= {a_n,a_n,c_n,c_n,a_n,b_n,c_n,b_n,b_n};

    const double x3D[27]= {a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
    const double y3D[27]= {a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
    const double z3D[27]= {a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};
    const double x3D_n[27]= {a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n};
    const double y3D_n[27]= {a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n};
    const double z3D_n[27]= {a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n};


//gaussian weights
// 1D --------------------------------
    _weight1[0][0]= 5./9.;
    _weight1[0][1]= 8./9.;
    _weight1[0][2]= 5./9.;

// 2D --------------------------------
    const double m=_weight1[0][0]*_weight1[0][1];
    const double l=_weight1[0][0]*_weight1[0][0];
    const double h=_weight1[0][1]*_weight1[0][1];
    const double weight[9]= {l,m,l,m,h,m,l,m,l};
    for ( int i=0; i<_NoGauss1[1]; i++ ) _weight1[1][i]= weight[i];

// 3D --------------------------------
    const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0];
    const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
    const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1];
    const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
    const double weight_3[27]= {w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
    for ( int i=0; i<_NoGauss1[2]; i++ ) _weight1[2][i]= weight_3[i];

#ifdef PRINT_TIME  //  TC +++++++++++++++ 
    std::clock_t start_time=std::clock();
#endif             //  TC +++++++++++++++ 

    /// General formula for test function: N_i = (1-0.5*|xi_i|)*((2|xi_i|-1)*xi*xi + xi_i*xi + (1-|xi_i|)) for xi, eta, chi
    /// xi is the gaussian point and xi_i is given by _CooE3[i] for the 1D case
    /// Loop i is over the shape functions. Loop j is over the gaussian points
    int nGauss, nShape, Dim;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _NoGauss1[0];
    nShape = _NoShape[0];
    Dim = 1;
    for ( int i=0; i<nShape; i++ ) {
        for ( int j=0; j<nGauss; j++ ) {
            double xx=x1D[j];
            double point[1];
            point[0] = xx;
            double xx_n=x1D_n[j];
            //shape functions
            _phi_map1[0][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );                 // phi
            _dphidxez_map1[0][i*nGauss+j] = Rec_Quad_DPhi ( i, point, Dim, 0 );        // dphi / dxi
            _dphidxx_map1[0][i*nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, 0, 0 );     // d2phi / d2eta
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _NoGauss1[1];
    nShape = _NoShape[1];
    Dim = 2;
    for ( int i=0; i<nShape; i++ ) {// loop over test function id
        for ( int j=0; j<nGauss; j++ ) {// loop over gauss nodes
            double point[2], pointNode[2];
            point[0] = x2D[j];
            point[1] = y2D[j];
            pointNode[0] = x2D_n[j];
            pointNode[1] = y2D_n[j];
            //shape functions
            _phi_map1[1][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );
            for ( int dir=0; dir<Dim; dir++ ) {
                _dphidxez_map1[1][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, point, Dim, dir );
                _dphidxez_map1_nodes[1][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, pointNode, Dim, dir );
                for ( int dir2=0; dir2<Dim; dir2++ )
                    _dphidxx_map1[1][ ( i+ ( dir*Dim + dir2 ) *nGauss ) *nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, dir, dir2 ); //d2/dxdy gaussian points
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nGauss = _NoGauss1[2];
        nShape = _NoShape[2];
	Dim = 3;
        // loop for test functions
        for ( int i=0; i<nShape; i++ ) {// loop over test function id
            for ( int j=0; j<nGauss; j++ ) {// loop over gauss nodes
                double point[3], pointNode[3];
                point[0] = x3D[j];
                point[1] = y3D[j];
                point[2] = z3D[j];
                pointNode[0] = x3D_n[j];
                pointNode[1] = y3D_n[j];
                pointNode[2] = z3D_n[j];
                //shape functions
                _phi_map1[2][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );
                for ( int dir=0; dir<Dim; dir++ ) {
                    _dphidxez_map1[2][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, point, Dim, dir );
                    _dphidxez_map1_nodes[2][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, pointNode, Dim, dir );
                    for ( int dir2=0; dir2<Dim; dir2++ )
                        _dphidxx_map1[2][ ( i+ ( dir*Dim + dir2 ) *nGauss ) *nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, dir, dir2 ); //d2/dxdy gaussian points
                }
            }
        }
    }
    return;
}

void MGFE::init_qua_Grec() {
//               ********************************************
//                               QUAD 9
// 		 ********************************************
//
// 			       3 ______6_____ 2
// 				|            |
// 			        |            |
//			       7|      8     |5
//				|            |
// 			        |____________|
//			       0       4      1

//gaussian coordinates
    const double a=-sqrt ( 3./5. );
    const double b=0.;
    const double c=-a;
    const double a_n=-1.;
    const double b_n=0.;
    const double c_n=1.;

    const double x1D[3]= {a,b,c};
    const double x1D_n[3]= {a_n,b_n,c_n};

    const double x2D[9]= {a,a,a,b,b,b,c,c,c};
    const double y2D[9]= {a,b,c,a,b,c,a,b,c};
    const double x2D_n[9]= {a_n,c_n,c_n,a_n,b_n,c_n,b_n,a_n,b_n};
    const double y2D_n[9]= {a_n,a_n,c_n,c_n,a_n,b_n,c_n,b_n,b_n};

    const double x3D[27]= {a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
    const double y3D[27]= {a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
    const double z3D[27]= {a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};
    const double x3D_n[27]= {a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,a_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,b_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n,c_n};
    const double y3D_n[27]= {a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n,a_n,a_n,a_n,b_n,b_n,b_n,c_n,c_n,c_n};
    const double z3D_n[27]= {a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n,a_n,b_n,c_n};


//gaussian weights
// 1D --------------------------------
    _Gweight1[3][0]= 5./9.;
    _Gweight1[3][1]= 8./9.;
    _Gweight1[3][2]= 5./9.;

// 2D --------------------------------
    const double m=_Gweight1[3][0]*_Gweight1[3][1];
    const double l=_Gweight1[3][0]*_Gweight1[3][0];
    const double h=_Gweight1[3][1]*_Gweight1[3][1];
    const double weight[9]= {l,m,l,m,h,m,l,m,l};
    for ( int i=0; i<_GNoGauss1[1][1]; i++ ) _Gweight1[4][i]= weight[i];

// 3D --------------------------------
    const double w1=_Gweight1[3][0]*_Gweight1[3][0]*_Gweight1[3][0];
    const double w2=_Gweight1[3][0]*_Gweight1[3][0]*_Gweight1[3][1];
    const double w3=_Gweight1[3][0]*_Gweight1[3][1]*_Gweight1[3][1];
    const double w4=_Gweight1[3][1]*_Gweight1[3][1]*_Gweight1[3][1];
    const double weight_3[27]= {w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
    for ( int i=0; i<_GNoGauss1[1][2]; i++ ) _Gweight1[5][i]= weight_3[i];

#ifdef PRINT_TIME  //  TC +++++++++++++++ 
    std::clock_t start_time=std::clock();
#endif             //  TC +++++++++++++++ 

    int nGauss, nShape, Dim;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _GNoGauss1[1][0];
    nShape = _GNoShape[1][0];
    Dim = 1;
    for ( int i=0; i<nShape; i++ ) {
        for ( int j=0; j<nGauss; j++ ) {
            double xx=x1D[j];
            double point[1];
            point[0] = xx;
            double xx_n=x1D_n[j];
            //shape functions
            _Gphi_map1[3][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );                 // phi
            _Gdphidxez_map1[3][i*nGauss+j] = Rec_Quad_DPhi ( i, point, Dim, 0 );        // dphi / dxi
            _Gdphidxx_map1[3][i*nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, 0, 0 );     // d2phi / d2eta
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _GNoGauss1[1][1];
    nShape = _GNoShape[1][1];
    Dim = 2;
    // loop for test functions
    for ( int i=0; i<nShape; i++ ) {// loop over test function id
        for ( int j=0; j<nGauss; j++ ) {// loop over gauss nodes
            double point[2], pointNode[2];
            point[0] = x2D[j];
            point[1] = y2D[j];
            pointNode[0] = x2D_n[j];
            pointNode[1] = y2D_n[j];
            //shape functions
            _Gphi_map1[4][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );
            for ( int dir=0; dir<2; dir++ ) {
                _Gdphidxez_map1[4][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, point, Dim, dir );
                _Gdphidxez_map1_nodes[4][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, pointNode, Dim, dir );
                for ( int dir2=0; dir2<2; dir2++ )
                    _Gdphidxx_map1[4][ ( i+ ( dir*2 + dir2 ) *nGauss ) *nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, dir, dir2 ); //d2/dxdy gaussian points
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nGauss = _GNoGauss1[1][2];
        nShape = _GNoShape[1][2];
	Dim = 3;
        // loop for test functions
        for ( int i=0; i<nShape; i++ ) {// loop over test function id
            for ( int j=0; j<nGauss; j++ ) {// loop over gauss nodes
                double point[3], pointNode[3];
                point[0] = x3D[j];
                point[1] = y3D[j];
                point[2] = z3D[j];
                pointNode[0] = x3D_n[j];
                pointNode[1] = y3D_n[j];
                pointNode[2] = z3D_n[j];
                //shape functions
                _Gphi_map1[5][i*nGauss+j] = Rec_Quad_Phi ( i, point, Dim );
                for ( int dir=0; dir<3; dir++ ) {
                    _Gdphidxez_map1[5][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, point, Dim, dir );
                    _Gdphidxez_map1_nodes[5][ ( i+dir*nGauss ) *nGauss + j] = Rec_Quad_DPhi ( i, pointNode, Dim, dir );
                    for ( int dir2=0; dir2<3; dir2++ )
                        _Gdphidxx_map1[5][ ( i+ ( dir*3 + dir2 ) *nGauss ) *nGauss+j] = Rec_Quad_D2Phi ( i, point, Dim, dir, dir2 ); //d2/dxdy gaussian points
                }
            }
        }
    }
    return;
}


void MGFE::init_qua_tri() {
    /*           ********************************************
    //                               TRI 6
    // 		 ********************************************
    //                              2
    // 				|\
    // 				| \
    // 				|  \
    // 			   r=l1 |   \ 1-s-r=l0
    //				|    \
    //				|     \
    // 			        |______\
    //			       0  s=l2 1
    */
//gaussian coordinates
// 1D --------------------------------
    const double x[3]= {-sqrt ( 3./5. ),0.,sqrt ( 3./5. ) };

// 2D --------------------------------
    const double alphabeta[4]= {1.-2.* ( 2./7. + sqrt ( 15. ) /21. ),  1.-2.* ( 2./7. - sqrt ( 15. ) /21. ), 2./7. + sqrt ( 15. ) /21., 2./7. - sqrt ( 15. ) /21.};
    const double r[7]= {1./3.,alphabeta[2],alphabeta[0],alphabeta[2],alphabeta[3],alphabeta[1],alphabeta[3]};
    const double s[7]= {1./3.,alphabeta[2],alphabeta[2],alphabeta[0],alphabeta[3],alphabeta[3],alphabeta[1]};

//3D --------------------------------
    const double p[14]= {0.31088591926330060980, 0.31088591926330060980,   1.-3.*0.31088591926330060980, 0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402,   1.-3.*0.092735250310891226402, 0.092735250310891226402,0.5-0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492};
    const double q[14]= {0.31088591926330060980, 1.-3.*0.31088591926330060980, 0.31088591926330060980,   0.31088591926330060980, 0.092735250310891226402, 1.-3.*0.092735250310891226402, 0.092735250310891226402,   0.092735250310891226402,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492};
    const double t[14]= {0.31088591926330060980, 0.31088591926330060980, 0.31088591926330060980,  1.-3.*0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402,  1.-3.*0.092735250310891226402,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.5-0.045503704125649649492};

//gaussian weights
// 1D --------------------------------
    _weight1[0][0]= 5./9.;
    _weight1[0][1]= 8./9.;
    _weight1[0][2]= 5./9.;

// 2D --------------------------------
    const double weight[7]= {9./80., 31./480. + sqrt ( 15. ) /2400., 31./480. + sqrt ( 15. ) /2400.,
                             31./480. + sqrt ( 15. ) /2400.,31./480. - sqrt ( 15. ) /2400., 31./480. - sqrt ( 15. ) /2400.,
                             31./480. - sqrt ( 15. ) /2400.
                            };

    for ( int i=0; i<_NoGauss1[1]; i++ ) _weight1[1][i]= weight[i];

// 3D --------------------------------
//  const double weight_3[5]={-16./120., 9./120., 9./120., 9./120., 9./120.};
    const double weight_3[14]= {0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730};
    for ( int i=0; i<_NoGauss1[2]; i++ ) _weight1[2][i]= weight_3[i];


// ****QUADRATIC SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates                  derivatives
//                                           |_______r_______|_______s_______
// phi0= 2.*(1-r-s)*((1.-r-s)-0.5)           | -3.+4.s+4.r   |  -3.+4.*s+4.*r
// phi1= 2.*r*(r-0.5)                        |    4.*r-1     |      0.
// phi2= 2.*s*(s-0.5)                        |    0.         |   4.*s-1
// phi3= 4.*r*(1.-r-s)                       |  4.-4.*s-8.*r |   -4.*r
// phi4= 4.*r*s                              |   4.*s        |   4.*r
// phi5= 4.*s*(1.-r-s)                       |  -4.*s        |   4-4.*r-8.*s
    int nShape, nGauss, Dim;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[0];
    nGauss = _NoGauss1[0];
    Dim = 1;
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[1];
            point[0] = x[i];
            const double xi=x[i];
            int off = nPhi*nGauss;
            _phi_map1[0][i + off]      = Rec_Quad_Phi ( nPhi, point, Dim );
            _dphidxez_map1[0][i + off] = Rec_Quad_DPhi ( nPhi, point, Dim, 0 );
            _dphidxx_map1[0][i + off]  = Rec_Quad_D2Phi ( nPhi, point, Dim, 0, 0 );
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[1];
    nGauss = _NoGauss1[1];
    for ( int nPhi=0; nPhi< nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = r[i];
            point[1] = s[i];
            _phi_map1[1][i+ ( nPhi ) *nGauss] = Tri_2d_QuadraticPhi ( nPhi, point );
            for ( int dir=0; dir<2; dir++ ) {
                int offset = dir*nShape + nPhi;
                _dphidxez_map1[1][i+ ( offset ) *nGauss] = Tri_2d_QuadraticDerPhi ( nPhi, point, dir );
                for ( int dir2=0; dir2<2; dir2++ ) {
                    _dphidxx_map1[1][i+ ( nPhi +6* ( dir*2 + dir2 ) ) *nGauss] = Tri_2d_QuadraticDer2Phi ( nPhi, point, dir, dir2 ); // d2 phi_i / d xi d xi
                }
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nShape = _NoShape[2];
        nGauss = _NoGauss1[2];
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = p[i];
                point[1] = q[i];
                point[2] = t[i];
                _phi_map1[2][i+ ( nPhi ) *nGauss] = Tri_3d_QuadraticPhi ( nPhi, point );
                for ( int dir1=0; dir1<3; dir1++ ) {
                    _dphidxez_map1[2][i+ ( nPhi + dir1*nShape ) *nGauss] = Tri_3d_QuadraticDerPhi ( nPhi, point, dir1 );
                    for ( int dir2=0; dir2<3; dir2++ )
                        _dphidxx_map1[2][i+ ( nPhi + ( dir2 + dir1*3 ) *nShape ) *nGauss] = Tri_3d_QuadraticDer2Phi ( nPhi, point, dir1, dir2 );
                }
            }
        }
    }
    return;
}

void MGFE::init_qua_Gtri() {
    /*           ********************************************
    //                               TRI 6
    // 		 ********************************************
    //                              2
    // 				|\
    // 				| \
    // 				|  \
    // 			   r=l1 |   \ 1-s-r=l0
    //				|    \
    //				|     \
    // 			        |______\
    //			       0  s=l2 1
    */
//gaussian coordinates
// 1D --------------------------------
    const double x[3]= {-sqrt ( 3./5. ),0.,sqrt ( 3./5. ) };

// 2D --------------------------------
    const double alphabeta[4]= {1.-2.* ( 2./7. + sqrt ( 15. ) /21. ),  1.-2.* ( 2./7. - sqrt ( 15. ) /21. ), 2./7. + sqrt ( 15. ) /21., 2./7. - sqrt ( 15. ) /21.};
    const double r[7]= {1./3.,alphabeta[2],alphabeta[0],alphabeta[2],alphabeta[3],alphabeta[1],alphabeta[3]};
    const double s[7]= {1./3.,alphabeta[2],alphabeta[2],alphabeta[0],alphabeta[3],alphabeta[3],alphabeta[1]};

//3D --------------------------------
    const double p[14]= {0.31088591926330060980, 0.31088591926330060980,   1.-3.*0.31088591926330060980, 0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402,   1.-3.*0.092735250310891226402, 0.092735250310891226402,0.5-0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492};
    const double q[14]= {0.31088591926330060980, 1.-3.*0.31088591926330060980, 0.31088591926330060980,   0.31088591926330060980, 0.092735250310891226402, 1.-3.*0.092735250310891226402, 0.092735250310891226402,   0.092735250310891226402,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492};
    const double t[14]= {0.31088591926330060980, 0.31088591926330060980, 0.31088591926330060980,  1.-3.*0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402,  1.-3.*0.092735250310891226402,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.5-0.045503704125649649492};

//gaussian weights
// 1D --------------------------------
    _Gweight1[0][0]= 5./9.;
    _Gweight1[0][1]= 8./9.;
    _Gweight1[0][2]= 5./9.;

// 2D --------------------------------
    const double weight[7]= {9./80., 31./480. + sqrt ( 15. ) /2400., 31./480. + sqrt ( 15. ) /2400.,
                             31./480. + sqrt ( 15. ) /2400.,31./480. - sqrt ( 15. ) /2400., 31./480. - sqrt ( 15. ) /2400.,
                             31./480. - sqrt ( 15. ) /2400.
                            };

    for ( int i=0; i<_GNoGauss1[0][1]; i++ ) _Gweight1[1][i]= weight[i];

// 3D --------------------------------
    const double weight_3[14]= {0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730};
    for ( int i=0; i<_GNoGauss1[0][2]; i++ ) _Gweight1[2][i]= weight_3[i];


// ****QUADRATIC SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates                  derivatives
//                                           |_______r_______|_______s_______
// phi0= 2.*(1-r-s)*((1.-r-s)-0.5)           | -3.+4.s+4.r   |  -3.+4.*s+4.*r
// phi1= 2.*r*(r-0.5)                        |    4.*r-1     |      0.
// phi2= 2.*s*(s-0.5)                        |    0.         |   4.*s-1
// phi3= 4.*r*(1.-r-s)                       |  4.-4.*s-8.*r |   -4.*r
// phi4= 4.*r*s                              |   4.*s        |   4.*r
// phi5= 4.*s*(1.-r-s)                       |  -4.*s        |   4-4.*r-8.*s
    int nShape, nGauss, Dim;
    Dim = 1;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nShape = _GNoShape[0][0];
    nGauss = _GNoGauss1[0][0];
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[1];
            point[0] = x[i];
            const double xi=x[i];
            int off = nPhi*nGauss;
            _Gphi_map1[0][i + off]      = Rec_Quad_Phi ( nPhi, point, Dim );
            _Gdphidxez_map1[0][i + off] = Rec_Quad_DPhi ( nPhi, point, Dim, 0 );
            _Gdphidxx_map1[0][i + off]  = Rec_Quad_D2Phi ( nPhi, point, Dim, 0, 0 );
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nShape = _GNoShape[0][1];
    nGauss = _GNoGauss1[0][1];
    for ( int nPhi=0; nPhi< nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = r[i];
            point[1] = s[i];
            _Gphi_map1[1][i+ ( nPhi ) *nGauss] = Tri_2d_QuadraticPhi ( nPhi, point );
            for ( int dir=0; dir<2; dir++ ) {
                int offset = dir*nShape + nPhi;
                _Gdphidxez_map1[1][i+ ( offset ) *nGauss] = Tri_2d_QuadraticDerPhi ( nPhi, point, dir );
                for ( int dir2=0; dir2<2; dir2++ ) {
                    _Gdphidxx_map1[1][i+ ( nPhi +6* ( dir*2 + dir2 ) ) *nGauss] = Tri_2d_QuadraticDer2Phi ( nPhi, point, dir, dir2 ); // d2 phi_i / d xi d xi
                }
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nShape = _GNoShape[0][2];
        nGauss = _GNoGauss1[0][2];
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = p[i];
                point[1] = q[i];
                point[2] = t[i];
                _Gphi_map1[2][i+ ( nPhi ) *nGauss] = Tri_3d_QuadraticPhi ( nPhi, point );
                for ( int dir1=0; dir1<3; dir1++ ) {
                    _Gdphidxez_map1[2][i+ ( nPhi + dir1*nShape ) *nGauss] = Tri_3d_QuadraticDerPhi ( nPhi, point, dir1 );
                    for ( int dir2=0; dir2<3; dir2++ )
                        _Gdphidxx_map1[2][i+ ( nPhi + ( dir2 + dir1*3 ) *nShape ) *nGauss] = Tri_3d_QuadraticDer2Phi ( nPhi, point, dir1, dir2 );
                }
            }
        }
    }
    return;
}


/// /// This function generates the Lagrangian linear shape functions
void MGFE::init_lin (
) {// ================================
    if ( _GlobalFE==1 ) {
        init_lin_Grec();
        init_lin_Gtri();
    }
#if ELTYPE==27
    init_lin_rec();
#endif

#if ELTYPE==10
    init_lin_tri();
#endif
    return;
}

void MGFE::init_lin_tri() {
    /* //        ********************************************
    //                               TRI 3
    // 		 ********************************************
    //                              2
    // 				|\
    // 				| \
    // 				|  \
    // 			   r=l1 |   \ 1-s-r=l0
    //				|    \
    //				|     \
    // 			        |______\
    //			       0  s=l2 1
    */
//gaussian coordinates
// 1D --------------------------------
    const double x[3]= {-sqrt ( 3./5. ),0.,sqrt ( 3./5. ) };
// 2D --------------------------------
    const double alphabeta[4]= {1.-2.* ( 2./7. + sqrt ( 15. ) /21. ),  1.-2.* ( 2./7. - sqrt ( 15. ) /21. ), 2./7. + sqrt ( 15. ) /21., 2./7. - sqrt ( 15. ) /21.};
    const double r[7]= {1./3.,alphabeta[2],alphabeta[0],alphabeta[2],alphabeta[3],alphabeta[1],alphabeta[3]};
    const double s[7]= {1./3.,alphabeta[2],alphabeta[2],alphabeta[0],alphabeta[3],alphabeta[3],alphabeta[1]};
//3D --------------------------------
    const double p[14]= {0.31088591926330060980, 0.31088591926330060980,   1.-3.*0.31088591926330060980, 0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402,   1.-3.*0.092735250310891226402, 0.092735250310891226402,0.5-0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492};
    const double q[14]= {0.31088591926330060980, 1.-3.*0.31088591926330060980, 0.31088591926330060980,   0.31088591926330060980, 0.092735250310891226402, 1.-3.*0.092735250310891226402, 0.092735250310891226402,   0.092735250310891226402,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492};
    const double t[14]= {0.31088591926330060980, 0.31088591926330060980, 0.31088591926330060980,  1.-3.*0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402,  1.-3.*0.092735250310891226402,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.5-0.045503704125649649492};

//gaussian weights
// 1D --------------------------------
    _weight1[0][0]= 5./9.;
    _weight1[0][1]= 8./9.;
    _weight1[0][2]= 5./9.;

// 2D --------------------------------
    const double weight[7]= {9./80., 31./480. + sqrt ( 15. ) /2400., 31./480. + sqrt ( 15. ) /2400.,
                             31./480. + sqrt ( 15. ) /2400.,31./480. - sqrt ( 15. ) /2400., 31./480. - sqrt ( 15. ) /2400.,
                             31./480. - sqrt ( 15. ) /2400.
                            };

    for ( int i=0; i<_NoGauss1[1]; i++ ) _weight1[1][i]= weight[i];


// 3D --------------------------------
    const double weight_3[14]= {0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730};
    for ( int i=0; i<_NoGauss1[2]; i++ ) _weight1[2][i]= weight_3[i];

// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates              derivatives
//                                      |_______r_______|_______s_______
// phi0= 1-r-s                          |      -1.      |      -1.
// phi1= r                              |       1.      |       0.
// phi2= s                              |       0.      |       1.
    int nShape, nGauss;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[0];
    nGauss = _NoGauss1[0];
    for ( int i=0; i<nGauss; i++ ) {
        const double xi=x[i];
        _phi_map1[0][i]=1.-xi;
        _phi_map1[0][i+nGauss]= 1.+ xi;
        _dphidxez_map1[0][i]=-1.;
        _dphidxez_map1[0][i+nGauss]=1.;
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[1];
    nGauss = _NoGauss1[1];
    for ( int nPhi=0; nPhi< nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = r[i];
            point[1] = s[i];
            _phi_map1[1][i+ ( nPhi ) *nGauss] = Tri_2d_LinearPhi ( nPhi, point );
            for ( int dir=0; dir<2; dir++ ) {
                int offset = dir*nShape + nPhi;
                _dphidxez_map1[1][i+ ( offset ) *nGauss] = Tri_2d_LinearDerPhi ( nPhi, point, dir );
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nShape = _NoShape[2];
        nGauss = _NoGauss1[2];
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = p[i];
                point[1] = q[i];
                point[2] = t[i];
                _phi_map1[2][i+ ( nPhi ) *nGauss] = Tri_3d_LinearPhi ( nPhi, point );
                for ( int dir1=0; dir1<3; dir1++ )
                    _dphidxez_map1[2][i+ ( nPhi + dir1*nShape ) *nGauss] = Tri_3d_LinearDerPhi ( nPhi, point, dir1 );
            }
        }
    }
    return;
}

void MGFE::init_lin_rec() {
//               ********************************************
//                               QUAD 4
// 		 ********************************************
//
// 			       3 ___________2
// 				|           |
// 			        |           |
//			        |           |
//				|           |
// 			        |___________|
//			       0             1

//gaussian coordinates
    const double a=-sqrt ( 3./5. );
    const double b=0.;
    const double c=-a;

    const double x1D[3]= {a,b,c};

    const double x2D[9]= {a,a,a,b,b,b,c,c,c};
    const double y2D[9]= {a,b,c,a,b,c,a,b,c};

    const double x3D[27]= {a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
    const double y3D[27]= {a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
    const double z3D[27]= {a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};


//gaussian weights
// 1D --------------------------------
    _weight1[0][0]= 5./9.;
    _weight1[0][1]= 8./9.;
    _weight1[0][2]= 5./9.;

// 2D --------------------------------
    const double m=_weight1[0][0]*_weight1[0][1];
    const double l=_weight1[0][0]*_weight1[0][0];
    const double h=_weight1[0][1]*_weight1[0][1];
    const double weight[9]= {l,m,l,m,h,m,l,m,l};
    for ( int i=0; i<_NoGauss1[1]; i++ ) _weight1[1][i]= weight[i];

// 3D --------------------------------
    const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0];
    const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
    const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1];
    const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
    const double weight_3[27]= {w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
    for ( int i=0; i<_NoGauss1[2]; i++ ) _weight1[2][i]= weight_3[i];


// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D                                            derivatives
//                                           |_______x_______|_______y_______
// phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// phi1= x*(1.-y)                            |     1.-y      |      -x
// phi2= x*y                                 |       y       |       x
// phi3= y*(1.-x)                            |      -y       |     1.-x

    int nShape, nGauss, Dim;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[0];
    nGauss = _NoGauss1[0];
    Dim = 1;
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[1];
            point[0] = x1D[i];
            _phi_map1[0][i+nPhi*nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                         // phi
            _dphidxez_map1[0][i+nPhi*nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, 0 );                // dphi / dxi
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nShape = _NoShape[1];
    nGauss = _NoGauss1[1];
    Dim = 2;
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = x2D[i];
            point[1] = y2D[i];
            _phi_map1[1][i+nPhi*nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                           // phi
            for ( int dir=0; dir<2; dir++ )
                _dphidxez_map1[1][i+ ( nPhi + dir*nShape ) *nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, dir );           // dphi / dxi
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nShape = _NoShape[2];
        nGauss = _NoGauss1[2];
	Dim = 3;
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = x3D[i];
                point[1] = y3D[i];
                point[2] = z3D[i];
                _phi_map1[2][i+ ( nPhi ) *nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                 // phi
                for ( int dir=0; dir<3; dir++ )
                    _dphidxez_map1[2][i+ ( nPhi + dir*nShape ) *nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, dir );
            }
        }
    }//end if 3D
    return;
}

void MGFE::init_lin_Gtri() {
    /* //        ********************************************
    //                               TRI 3
    // 		 ********************************************
    //                          2
    // 				|\
    // 				| \
    // 				|  \
    // 			   r=l1 |   \ 1-s-r=l0
    //				|    \
    //				|     \
    // 			        |______\
    //			       0  s=l2  1
    */
//gaussian coordinates
// 1D --------------------------------
    const double x[3]= {-sqrt ( 3./5. ),0.,sqrt ( 3./5. ) };

// 2D --------------------------------
    const double alphabeta[4]= {1.-2.* ( 2./7. + sqrt ( 15. ) /21. ),  1.-2.* ( 2./7. - sqrt ( 15. ) /21. ), 2./7. + sqrt ( 15. ) /21., 2./7. - sqrt ( 15. ) /21.};
    const double r[7]= {1./3.,alphabeta[2],alphabeta[0],alphabeta[2],alphabeta[3],alphabeta[1],alphabeta[3]};
    const double s[7]= {1./3.,alphabeta[2],alphabeta[2],alphabeta[0],alphabeta[3],alphabeta[3],alphabeta[1]};

//3D --------------------------------
    const double p[14]= {0.31088591926330060980, 0.31088591926330060980,   1.-3.*0.31088591926330060980, 0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402,   1.-3.*0.092735250310891226402, 0.092735250310891226402,0.5-0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492};
    const double q[14]= {0.31088591926330060980, 1.-3.*0.31088591926330060980, 0.31088591926330060980,   0.31088591926330060980, 0.092735250310891226402, 1.-3.*0.092735250310891226402, 0.092735250310891226402,   0.092735250310891226402,0.5-0.045503704125649649492,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492};
    const double t[14]= {0.31088591926330060980, 0.31088591926330060980, 0.31088591926330060980,  1.-3.*0.31088591926330060980, 0.092735250310891226402, 0.092735250310891226402, 0.092735250310891226402,  1.-3.*0.092735250310891226402,0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.045503704125649649492,0.5-0.045503704125649649492,0.5-0.045503704125649649492};

//gaussian weights
// 1D --------------------------------
    _Gweight1[0][0]= 5./9.;
    _Gweight1[0][1]= 8./9.;
    _Gweight1[0][2]= 5./9.;

// 2D --------------------------------
    const double weight[7]= {9./80., 31./480. + sqrt ( 15. ) /2400., 31./480. + sqrt ( 15. ) /2400.,
                             31./480. + sqrt ( 15. ) /2400.,31./480. - sqrt ( 15. ) /2400., 31./480. - sqrt ( 15. ) /2400.,
                             31./480. - sqrt ( 15. ) /2400.
                            };

    for ( int i=0; i<_GNoGauss1[0][1]; i++ ) _Gweight1[1][i]= weight[i];


// 3D --------------------------------
    const double weight_3[14]= {0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.018781320953002641800,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.012248840519393658257,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730,0.0070910034628469110730};
    for ( int i=0; i<_GNoGauss1[0][2]; i++ ) _Gweight1[2][i]= weight_3[i];

// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D in triangular coordinates              derivatives
//                                      |_______r_______|_______s_______
// phi0= 1-r-s                          |      -1.      |      -1.
// phi1= r                              |       1.      |       0.
// phi2= s                              |       0.      |       1.
    int nShape, nGauss;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nShape = _GNoShape[0][0];
    nGauss = _GNoGauss1[0][0];
    for ( int i=0; i<nGauss; i++ ) {
        const double xi=x[i];
        _Gphi_map1[0][i]=1.-xi;
        _Gphi_map1[0][i+nGauss]= 1.+ xi;
        _Gdphidxez_map1[0][i]=-1.;
        _Gdphidxez_map1[0][i+nGauss]=1.;
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nShape = _GNoShape[0][1];
    nGauss = _GNoGauss1[0][1];
    for ( int nPhi=0; nPhi< nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = r[i];
            point[1] = s[i];
            _Gphi_map1[1][i+ ( nPhi ) *nGauss] = Tri_2d_LinearPhi ( nPhi, point );
            for ( int dir=0; dir<2; dir++ ) {
                int offset = dir*nShape + nPhi;
                _Gdphidxez_map1[1][i+ ( offset ) *nGauss] = Tri_2d_LinearDerPhi ( nPhi, point, dir );
            }
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nShape = _GNoShape[0][2];
        nGauss = _GNoGauss1[0][2];
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = p[i];
                point[1] = q[i];
                point[2] = t[i];
                _Gphi_map1[2][i+ ( nPhi ) *nGauss] = Tri_3d_LinearPhi ( nPhi, point );
                for ( int dir1=0; dir1<3; dir1++ ) _Gdphidxez_map1[2][i+ ( nPhi + dir1*nShape ) *nGauss] = Tri_3d_LinearDerPhi ( nPhi, point, dir1 );
            }
        }
    }

    return;
}


void MGFE::init_lin_Grec() {
//               ********************************************
//                               QUAD 4
// 		 ********************************************
//
// 			       3 ___________2
// 				|           |
// 			        |           |
//			        |           |
//				|           |
// 			        |___________|
//			       0             1

//gaussian coordinates
    const double a=-sqrt ( 3./5. );
    const double b=0.;
    const double c=-a;

    const double x1D[3]= {a,b,c};

    const double x2D[9]= {a,a,a,b,b,b,c,c,c};
    const double y2D[9]= {a,b,c,a,b,c,a,b,c};

    const double x3D[27]= {a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
    const double y3D[27]= {a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
    const double z3D[27]= {a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};


//gaussian weights
// 1D --------------------------------
    _Gweight1[3][0]= 5./9.;
    _Gweight1[3][1]= 8./9.;
    _Gweight1[3][2]= 5./9.;

// 2D --------------------------------
    const double m=_Gweight1[3][0]*_Gweight1[3][1];
    const double l=_Gweight1[3][0]*_Gweight1[3][0];
    const double h=_Gweight1[3][1]*_Gweight1[3][1];
    const double weight[9]= {l,m,l,m,h,m,l,m,l};
    for ( int i=0; i<_GNoGauss1[1][1]; i++ ) _Gweight1[4][i]= weight[i];

// 3D --------------------------------
    const double w1=_Gweight1[3][0]*_Gweight1[3][0]*_Gweight1[3][0];
    const double w2=_Gweight1[3][0]*_Gweight1[3][0]*_Gweight1[3][1];
    const double w3=_Gweight1[3][0]*_Gweight1[3][1]*_Gweight1[3][1];
    const double w4=_Gweight1[3][1]*_Gweight1[3][1]*_Gweight1[3][1];
    const double weight_3[27]= {w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
    for ( int i=0; i<_GNoGauss1[1][2]; i++ ) _Gweight1[5][i]= weight_3[i];


// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D                                            derivatives
//                                           |_______x_______|_______y_______
// phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// phi1= x*(1.-y)                            |     1.-y      |      -x
// phi2= x*y                                 |       y       |       x
// phi3= y*(1.-x)                            |      -y       |     1.-x

    int nShape, nGauss, Dim;
// -----------------------------------------------------------------------
//                            1D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _GNoGauss1[1][0];
    nShape = _GNoShape[1][0];
    Dim = 1;
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[1];
            point[0] = x1D[i];
            _Gphi_map1[3][i+nPhi*nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                         // phi
            _Gdphidxez_map1[3][i+nPhi*nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, 0 );                // dphi / dxi
        }
    }
// -----------------------------------------------------------------------
//                            2D ELEMENT
// -----------------------------------------------------------------------
    nGauss = _GNoGauss1[1][1];
    nShape = _GNoShape[1][1];
    Dim = 2;
    for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
        for ( int i=0; i<nGauss; i++ ) {
            double point[2];
            point[0] = x2D[i];
            point[1] = y2D[i];
            _Gphi_map1[4][i+nPhi*nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                           // phi
            for ( int dir=0; dir<2; dir++ )
                _Gdphidxez_map1[4][i+ ( nPhi + dir*nShape ) *nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, dir );          // dphi / dxi
        }
    }
// -----------------------------------------------------------------------
//                            3D ELEMENT
// -----------------------------------------------------------------------
    if ( _dim == 3 ) {
        nGauss = _GNoGauss1[1][2];
        nShape = _GNoShape[1][2];
	Dim = 3;
        for ( int nPhi=0; nPhi<nShape; nPhi++ ) {
            for ( int i=0; i<nGauss; i++ ) {
                double point[3];
                point[0] = x3D[i];
                point[1] = y3D[i];
                point[2] = z3D[i];
                _Gphi_map1[5][i+ ( nPhi ) *nGauss] = Rec_Lin_Phi ( nPhi, point, Dim );                 // phi
                for ( int dir=0; dir<3; dir++ )
                    _Gdphidxez_map1[5][i+ ( nPhi + dir*nShape ) *nGauss] = Rec_Lin_DPhi ( nPhi, point, Dim, dir );          // dphi / dxi
            }
        }
    }//end if 3D
    return;
}

/// /// This function generates the Lagrangian piecewise shape functions
void MGFE::init_pie (
) {// ================================

//   CONTROLLARE PER IMPLEMENTAZIONE ELEMENTI MISTI
//
//
//
//               ********************************************
//                               S_0
//               ********************************************
//
//                               ___________
//                              |           |
//                              |     0     |
//                              |     x     |
//                              |           |
//                              |___________|
//

//gaussian coordinates
    const double a=-sqrt ( 3./5. );
    const double b=0.;
    const double c=-a;

    const double x1D[3]= {a,b,c};

    const double x2D[9]= {a,a,a,b,b,b,c,c,c};
    const double y2D[9]= {a,b,c,a,b,c,a,b,c};

    const double x3D[27]= {a,a,a,a,a,a,a,a,a,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,c,c};
    const double y3D[27]= {a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c,a,a,a,b,b,b,c,c,c};
    const double z3D[27]= {a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c,a,b,c};


//gaussian weights
// 1D --------------------------------
    _weight1[0][0]= 5./9.;
    _weight1[0][1]= 8./9.;
    _weight1[0][2]= 5./9.;

// 2D --------------------------------
    const double m=_weight1[0][0]*_weight1[0][1];
    const double l=_weight1[0][0]*_weight1[0][0];
    const double h=_weight1[0][1]*_weight1[0][1];
    const double weight[9]= {l,m,l,m,h,m,l,m,l};
    for ( int i=0; i<_NoGauss1[1]; i++ ) _weight1[1][i]= weight[i];

// 3D --------------------------------
    const double w1=_weight1[0][0]*_weight1[0][0]*_weight1[0][0];
    const double w2=_weight1[0][0]*_weight1[0][0]*_weight1[0][1];
    const double w3=_weight1[0][0]*_weight1[0][1]*_weight1[0][1];
    const double w4=_weight1[0][1]*_weight1[0][1]*_weight1[0][1];
    const double weight_3[27]= {w1,w2,w1,w2,w3,w2,w1,w2,w1,w2,w3,w2,w3,w4,w3,w2,w3,w2,w1,w2,w1,w2,w3,w2,w1,w2,w1};
    for ( int i=0; i<_NoGauss1[2]; i++ ) _weight1[2][i]= weight_3[i];


// ****LINEAR SHAPES AND DERIVATIVES****
// shape 2D                                            derivatives
//                                           |_______x_______|_______y_______
// phi0= (1.-x)*(1.-y)                       |     y.-1      |     x.-1
// phi1= x*(1.-y)                            |     1.-y      |      -x
// phi2= x*y                                 |       y       |       x
// phi3= y*(1.-x)                            |      -y       |     1.-x

//1D --------------------------------
    for ( int i=0; i<_NoGauss1[0]; i++ ) {
        _phi_map1[0][i]=1.;
        _dphidxez_map1[0][i]=0.;

        if ( NDOF_K>1 ) {
            _phi_map1[0][i+_NoGauss1[0]]= x1D[i];
            _dphidxez_map1[0][i+_NoGauss1[0]]=1.;
        }
    }

    // 2D --------------------------------
    for ( int i=0; i<_NoGauss1[1]; i++ ) {
        const double xx=x2D[i];
        const double yy=y2D[i];
        //shape functions
        _phi_map1[1][i+0*_NoGauss1[1]]=1.;
        // derivatives
        _dphidxez_map1[1][i+ ( 0 ) *_NoGauss1[1]]= 0.; //d/dx
        _dphidxez_map1[1][i+ ( 1 ) *_NoGauss1[1]]= 0.; //d/dy

        if ( NDOF_K>1 ) {
            _phi_map1[1][i+0*_NoGauss1[1]]=1.;
            _phi_map1[1][i+1*_NoGauss1[1]]=xx;
            _phi_map1[1][i+2*_NoGauss1[1]]=yy;

            //d/dx
            _dphidxez_map1[1][i+ ( 0 ) *_NoGauss1[1]]= 0.;
            _dphidxez_map1[1][i+ ( 1 ) *_NoGauss1[1]]= 1.;
            _dphidxez_map1[1][i+ ( 2 ) *_NoGauss1[1]]= 0.;
            //d/dy
            _dphidxez_map1[1][i+ ( 3 ) *_NoGauss1[1]]= 0.;
            _dphidxez_map1[1][i+ ( 4 ) *_NoGauss1[1]]= 0.;
            _dphidxez_map1[1][i+ ( 5 ) *_NoGauss1[1]]= 1.;

        }

    }

    // 3D -----------------------------------------------
    if ( _dim == 3 ) {

        for ( int i=0; i<_NoGauss1[2]; i++ ) {
            const double xx=x3D[i];
            const double yy=y3D[i];
            const double zz=z3D[i];

            //shape functions
            _phi_map1[2][i+ ( 0 ) *_NoGauss1[2]]= 1.;
            // derivatives
            _dphidxez_map1[2][i+ ( 0 ) *_NoGauss1[2]]= 0.; //d/dx
            _dphidxez_map1[2][i+ ( 1 ) *_NoGauss1[2]]= 0.; //d/dy
            _dphidxez_map1[2][i+ ( 2 ) *_NoGauss1[2]]= 0.; //d/dz


            if ( NDOF_K>1 ) {
                //shape functions
                _phi_map1[2][i+ ( 0 ) *_NoGauss1[2]]= 1.;
                _phi_map1[2][i+ ( 1 ) *_NoGauss1[2]]= xx;
                _phi_map1[2][i+ ( 2 ) *_NoGauss1[2]]= yy;
                _phi_map1[2][i+ ( 3 ) *_NoGauss1[2]]= zz;

                //d/dx
                _dphidxez_map1[2][i+ ( 0 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 1 ) *_NoGauss1[2]]= 1.;
                _dphidxez_map1[2][i+ ( 2 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 3 ) *_NoGauss1[2]]= 0.;
                //d/dy
                _dphidxez_map1[2][i+ ( 4 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 5 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 6 ) *_NoGauss1[2]]= 1.;
                _dphidxez_map1[2][i+ ( 7 ) *_NoGauss1[2]]= 0.;
                //d/dz
                _dphidxez_map1[2][i+ ( 8 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 9 ) *_NoGauss1[2]]= 0.;
                _dphidxez_map1[2][i+ ( 10 ) *_NoGauss1[2]]=0.;
                _dphidxez_map1[2][i+ ( 11 ) *_NoGauss1[2]]= 1.;
            }
        }
    }//end if 3D
    return;
}




// ================================
/// This function writes FEM
void MGFE::write (
    const std::string& name,// file <-
    const int kdim         // dimension <-
) {// ================================
    std::ofstream in ( name.c_str() );
    this->write_c ( in,kdim );
}

// ====================================================
/// This function writes shape and derivative values at the gaussian points
void MGFE::write_c (
    std::ostream& out,    // file <-
    const int kdim       // dimension  <-
) { // ================================================

    if ( !out ) {
        std::cout<<" Gauss Outfile "<<out<<" not opened."<<std::endl;
        exit ( 3 );
    }

    // heading
    out<< "gpoints" <<  _NoGauss1[kdim-1]  << std::endl  ;
    // max_nodes (decreasing level order)
    for ( int k=0; k<_NoGauss1[kdim-1]; k++ ) {
        out<< " weight " << std::setprecision ( 20 ) <<  _weight1[kdim-1][k] << std::endl;
        for ( int s=0; s<_NoShape[kdim-1]; s++ ) {
            out<<std::setprecision ( 20 ) <<_phi_map1[kdim-1][s*_NoGauss1[kdim-1]+k]<<" ";
            for ( int idim=0; idim<3; idim++ ) {
                out<<std::setprecision ( 20 ) <<
                _dphidxez_map1[kdim-1][ ( idim*_NoShape[kdim-1]+s ) *_NoGauss1[kdim-1]+k] << "  ";
            }
            out<< "  \n";
        }
        out<< "  \n";
    }
#ifdef PRINT_INFO
    std::cout << " MGFE::write_c:  dim " << kdim<< " fem with " << _NoGauss1[kdim-1] << " gaussian points and "
    << _NoShape[kdim-1] << " shape functions \n" << std::endl;
#endif
    return;
}


// ===========================================================
/// This function read the gaussian points
void MGFE::read_c (
    std::istream& infile, // file <-
    const int kdim       // dimension <-
) {// =========================================================


    // file
    if ( !infile ) {
        std::cout << " Read_c3:Gauss Input file not opened \n";
        std::exit ( 3 );
    }

    const int  bufLen = 30;
    char  buf[bufLen+1];
    sprintf ( buf,"0" );
    int ng;
    double dummy;
    int ngss=_NoGauss1[kdim-1];
    int nsh=_NoShape[kdim-1];

    // Reading file
    while ( strncmp ( buf,"gpoints",7 ) != 0 )    infile >> buf;
    // # of gaussian points
    infile >> ng;
    if ( _NoGauss1[kdim-1] != ng ) std::cout <<
        " MGFE::read_c: error _NoGauss is " << _NoGauss1[kdim-1] << std::endl;

    // Reading 3d shape functions at guassian points

    for ( int k=0; k<ngss; k++ ) {         // gauss points
        infile >> buf;
        infile >> dummy;
        _weight1[kdim-1][k]=dummy;
//         std::cout<<" "<< kdim-1 << " " << k << " " << dummy<<"\n";
        for ( int s=0; s< nsh; s++ ) {       // shape functions
            infile >> dummy;
            _phi_map1[kdim-1][s*ngss+k]=dummy;

            for ( int idim=0; idim<kdim; idim++ ) { // derivatives
                infile >> dummy;
                _dphidxez_map1[kdim-1][ ( idim*nsh+s ) *ngss+k]=dummy;
            }
        }
    }
    std::cout<<"\n\n**************************/*phi*/*******************************\n\n";
    for ( int k=0; k<ngss*nsh; k++ ) {
//   if (k%27==0) std::cout<< "phi" << k/27  << "\n\n";
        std::cout/*<< k<< " "*/ << _phi_map1[2][k] << " \n ";
//     std::cout <<  _weight1[2][k] << "\n";
    }

//   std::cout<<"\n\n**************************/*DERIVATA*/*******************************\n\n";
//   for (int i=0;i<ngss*nsh*3;i++) {
//        if (i%27==0) std::cout<< "dphi" << i/27  << "\n\n";
//     std::cout<<" "<< i<< " "<< " " <<_dphidxez_map1[2][i]<<"\n";
//   }

#ifdef PRINT_INFO
    std::cout << " MGFE::read_c:  dim " << kdim<< " fem with " << ngss << " gaussian points and "
    << nsh << " shape functions " << std::endl;
#endif

    return;
}

// ========================================
/// This function computes the  derivatives at
/// the gaussian point ng:
double MGFE::Jac ( //  jacobean ->
    const int ng, // gaussian point <-
    double x[],     // coordinates  <-
    double InvJac[] // Inverted Jacobean ->
) {// =====================================
    double Jac=0.;
#if  DIMENSION==1
    Jac = Jac1D ( ng,x,InvJac );
#endif
#if  DIMENSION==2
    Jac = Jac2D ( ng,x,InvJac );
#endif
#if  DIMENSION==3
    Jac = Jac3D ( ng,x,InvJac );
#endif
    return Jac;
}

// ========================================
/// This function computes the 3D derivatives at
/// the gaussian point ng:
double MGFE::Jac3D ( // 3D  jacobean ->
    const int ng,   // gaussian point <-
    double xyz[],     // coordinates <-
    double InvJac[] // Jacobean ->
) { // =====================================

    double x_xi=0., x_eta =0.,x_zeta =0.;
    double y_xi =0.,y_eta=0., y_zeta =0.;
    double z_xi=0., z_eta =0.,z_zeta =0.;
    int nshape=_NoShape[2];
    int offset=nshape*_NoGauss1[2];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[2]+ng;
        double dphidx=_dphidxez_map1[2][sng];
        double dphideta=_dphidxez_map1[2][sng+offset];
        double dphidzeta=_dphidxez_map1[2][sng+2*offset];

        x_xi  += xyz[s]*dphidx;
        x_eta += xyz[s]*dphideta;
        x_zeta +=xyz[s]*dphidzeta;
        y_xi  += xyz[s+nshape]*dphidx;
        y_eta += xyz[s+nshape]*dphideta;
        y_zeta +=xyz[s+nshape]*dphidzeta;
        z_xi  += xyz[s+2*nshape]*dphidx;
        z_eta += xyz[s+2*nshape]*dphideta;
        z_zeta +=xyz[s+2*nshape]*dphidzeta;

    }
    double det=x_xi* ( y_eta*z_zeta-y_zeta*z_eta )-
               x_eta* ( y_xi*z_zeta-y_zeta*z_xi ) +
               x_zeta* ( y_xi*z_eta-y_eta*z_xi );
    double invdet=1./det;
    InvJac[0]= ( y_eta*z_zeta-y_zeta*z_eta ) *invdet;
    InvJac[3]=- ( x_eta*z_zeta-x_zeta*z_eta ) *invdet; //
    InvJac[6]= ( y_zeta*x_eta-x_zeta*y_eta ) *invdet;
    InvJac[1]=- ( y_xi*z_zeta-y_zeta*z_xi ) *invdet; //
    InvJac[4]= ( x_xi*z_zeta-x_zeta*z_xi ) *invdet;
    InvJac[7]=- ( y_zeta*x_xi-y_xi*x_zeta ) *invdet; //
    InvJac[2]= ( z_eta*y_xi-y_eta*z_xi ) *invdet;
    InvJac[5]=- ( x_xi*z_eta-x_eta*z_xi ) *invdet; //
    InvJac[8]= ( x_xi*y_eta-y_xi*x_eta ) *invdet;


    return ( det );
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac1D ( // 2D  jacobean ->
    const int ng,   // Gaussian point <-
    double x[],       // coordinates  <-
    double InvJac[] // Jacobean
)  {// ===================================================

    double x_xi=0.;
    int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss1[0];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[0]+ng;
        x_xi  +=x[s]*_dphidxez_map1[0][sng];
//     x_eta +=x[s]*_dphidxez_map1[1][sng+offset];
//     y_xi  +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng];
//     y_eta +=x[s+NDOF_FEM]*_dphidxez_map1[1][sng+offset];
    }
    double det=x_xi;
    double idet=1./det;
//   InvJac[0]=y_eta*idet;    // dxi dx
//   InvJac[1]=-y_xi*idet;    //deta dx
//   InvJac[2]=-x_eta*idet;   // dxi dy
    InvJac[0]=x_xi*idet;     //deta dy

    return ( det );
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac2D ( // 2D  jacobean ->
    const int ng,   // Gaussian point <-
    double x[],       // coordinates  <-
    double InvJac[] // Inverted Jacobean
)  {// ===================================================

    double x_xi=0., x_eta =0.,y_xi =0., y_eta=0.;
    int nshape=_NoShape[1];
    int offset=_NoShape[1]*_NoGauss1[1];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[1]+ng;
        x_xi  +=x[s]*_dphidxez_map1[1][sng];
        x_eta +=x[s]*_dphidxez_map1[1][sng+offset];
        y_xi  +=x[s+nshape]*_dphidxez_map1[1][sng];
        y_eta +=x[s+nshape]*_dphidxez_map1[1][sng+offset];
    }
    double det=x_xi*y_eta-y_xi*x_eta;
    double idet=1./det;
    InvJac[0]=y_eta*idet;    // dxi dx
    InvJac[1]=-y_xi*idet;    //deta dx
    InvJac[2]=-x_eta*idet;   // dxi dy
    InvJac[3]=x_xi*idet;     //deta dy

    return ( det );
}

// ========================================
/// This function computes the  derivatives at
/// the gaussian point ng:
double MGFE::JacG ( //  jacobean ->
    const int ng, // gaussian point <-
    double x[],     // coordinates  <-
    double InvJac[], // Inverted Jacobean ->
    int FamilyType,
    int dim
) {// =====================================
    double Jac=0.;
    int nShape=_GNoShape[FamilyType][dim-1];
    int offset=_GNoShape[FamilyType][dim-1]*_GNoGauss1[FamilyType][dim-1];
    double *CoordsDer = new double[dim*dim];
    for ( int i=0; i<dim*dim; i++ ) CoordsDer[i] = 0.;
    for ( int s=0; s<nShape; s++ ) {
        int sng=s*_GNoGauss1[FamilyType][dim-1]+ng;
        for ( int row=0; row<dim; row++ ) {
            for ( int col=0; col<dim; col++ ) {
                CoordsDer[row*dim + col] += x[s + row*nShape] * _Gdphidxez_map1[dim-1 + 3*FamilyType][sng+col*offset];
            }
        }
    }
    double det = ComputeInverseMatrix(CoordsDer, InvJac, dim);
    delete [] CoordsDer;
    return ( det );
    return Jac;
}

void MGFE::JacOnGivenCanCoords (
    int dim,
    double ElemCoords[],
    double CanCoords[],
    double InvJac[],
    int FamilyType,
    int nShape
) {
    double *LocDPhi = new double[dim];
    double *CoordsDer = new double[dim*dim];

    for ( int dir=0; dir<dim*dim; dir++ ) CoordsDer[dir] = 0.;
    for ( int s=0; s<nShape; s++ ) {
        for ( int dir=0; dir<dim; dir++ ) LocDPhi[dir] = FirstDerivateOfLocalPhi ( s, CanCoords, dim, dir, FamilyType );

        for ( int dir1=0; dir1<dim; dir1++ )
            for ( int dir2=0; dir2<dim; dir2++ )
                CoordsDer[dir1*dim + dir2] += ElemCoords[s + dir1*nShape] * LocDPhi[dir2];
    }

    double det=0.;
    if ( dim==3 ) {
        for ( int dir=0; dir<dim; dir++ ) {
            int sign = 1 - 2* ( dir%2 );
            int idx1 = ( dir +1 ) %dim;
            int idx2 = ( dir +2 ) %dim;
            det += sign* ( CoordsDer[idx1*dim + idx1]*CoordsDer[idx2*dim + idx2] - CoordsDer[idx1*dim + idx2]*CoordsDer[idx2*dim + idx1] );
        }
    }
    if ( dim==2 ) {
        det = ( CoordsDer[0]*CoordsDer[3] - CoordsDer[1]*CoordsDer[2] );
    }
    double idet=1./det;

    if ( dim==2 ) {
        InvJac[0] =  CoordsDer[3]*idet;    // dxi dx
        InvJac[1] = -CoordsDer[2]*idet;    // deta dx
        InvJac[2] = -CoordsDer[1]*idet;    // dxi dy
        InvJac[3] =  CoordsDer[0]*idet;    // deta dy
    }
    if ( dim==3 ) {
        for ( int row=0; row<dim; row++ ) {// LOOP OVER ROWS - row1<row2
            int sign_row =  1 - 2* ( row%2 );
            int row1 = ( ( row +1 ) %dim < ( row +2 ) %dim ) ? ( row +1 ) %dim : ( row +2 ) %dim;
            int row2 = ( ( row +1 ) %dim > ( row +2 ) %dim ) ? ( row +1 ) %dim : ( row +2 ) %dim;
            for ( int col=0; col<dim; col++ ) {// LOOP OVER COLUMNS  - col1<col2
                int col1 = ( ( col +1 ) %dim < ( col +2 ) %dim ) ? ( col +1 ) %dim : ( col +2 ) %dim;
                int col2 = ( ( col +1 ) %dim > ( col +2 ) %dim ) ? ( col +1 ) %dim : ( col +2 ) %dim;
                int sign_col =  1 - 2* ( col%2 );
                InvJac[row*dim + col] = sign_col*sign_row * idet *
                                        ( CoordsDer[row1*dim + col1]*CoordsDer[row2*dim + col2] - CoordsDer[row1*dim + col2]*CoordsDer[row2*dim + col1] );
            }
        }
    }

    delete [] LocDPhi;
    delete [] CoordsDer;
    return;
}


void MGFE::get_dphi_on_given_node (
    const int dim,
    double ElemCoords[],
    double CanPos[],
    double dphi[]
) {
    double *InvJac    = new double[dim*dim];
    double *gradphi_g = new double[dim];
    int nShape=_NoShape[dim-1];
    JacOnGivenCanCoords ( dim, ElemCoords, CanPos, InvJac, _FamilyType, nShape );
    for ( int eln=0; eln<nShape; eln++ )    {
        for ( int idim=0; idim<dim; idim++ )
            gradphi_g[idim] = FirstDerivateOfLocalPhi ( eln, CanPos, dim, idim, _FamilyType );

        for ( int idim=0; idim<dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<dim; jdim++ ) sum += InvJac[jdim+idim*dim]*gradphi_g[jdim];
            dphi[eln+idim* nShape] = sum;
        }
    }

    delete [] gradphi_g;
    delete [] InvJac;
    return;
}

void MGFE::get_dphi_on_given_nodeG (
    const int dim,
    double ElemCoords[],
    double CanPos[],
    double dphi[],
    int FamilyType
) {
    double *InvJac    = new double[dim*dim];
    double *gradphi_g = new double[dim];
    int nShape = _GNoShape[FamilyType][dim-1];
    JacOnGivenCanCoords ( dim, ElemCoords, CanPos, InvJac, FamilyType, nShape );

    for ( int eln=0; eln<nShape; eln++ )    {
        for ( int idim=0; idim<dim; idim++ )
            gradphi_g[idim] = FirstDerivateOfLocalPhi ( eln, CanPos, dim, idim, FamilyType );

        for ( int idim=0; idim<dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<dim; jdim++ ) sum += InvJac[jdim+idim*dim]*gradphi_g[jdim];
            dphi[eln+idim* nShape] = sum;
        }
    }

    delete [] gradphi_g;
    delete [] InvJac;
    return;
}
// ========================================
/// This function computes the  derivatives at
/// the nodal point:
double MGFE::Jac_nodes ( //  jacobean ->
    const int ng, // nodal point <-
    double x[],     // coordinates  <-
    double InvJac[] // Jacobean ->
) {// =====================================
    double Jac=0.;
#if  DIMENSION==1
    Jac = Jac1D_nodes ( ng,x,InvJac );
#endif
#if  DIMENSION==2
    Jac = Jac2D_nodes ( ng,x,InvJac );
#endif
#if  DIMENSION==3
    Jac = Jac3D_nodes ( ng,x,InvJac );
#endif
    return Jac;
}

// ========================================
/// This function computes the 3D derivatives at
/// the nodal point ng:
double MGFE::Jac3D_nodes ( // 3D  jacobean ->
    const int ng,   // nodal point <-
    double xyz[],     // coordinates <-
    double InvJac[] // Jacobean ->
) { // =====================================

    double x_xi=0., x_eta =0.,x_zeta =0.;
    double y_xi =0.,y_eta=0., y_zeta =0.;
    double z_xi=0., z_eta =0.,z_zeta =0.;
    int nshape=_NoShape[2];
    int offset=nshape*_NoGauss1[2];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[2]+ng;
        double dphidx=_dphidxez_map1_nodes[2][sng];
        double dphideta=_dphidxez_map1_nodes[2][sng+offset];
        double dphidzeta=_dphidxez_map1_nodes[2][sng+2*offset];

        x_xi  += xyz[s]*dphidx;
        x_eta += xyz[s]*dphideta;
        x_zeta +=xyz[s]*dphidzeta;
        y_xi  += xyz[s+nshape]*dphidx;
        y_eta += xyz[s+nshape]*dphideta;
        y_zeta +=xyz[s+nshape]*dphidzeta;
        z_xi  += xyz[s+2*nshape]*dphidx;
        z_eta += xyz[s+2*nshape]*dphideta;
        z_zeta +=xyz[s+2*nshape]*dphidzeta;

    }
    double det=x_xi* ( y_eta*z_zeta-y_zeta*z_eta )-
               x_eta* ( y_xi*z_zeta-y_zeta*z_xi ) +
               x_zeta* ( y_xi*z_eta-y_eta*z_xi );
    double invdet=1./det;
    InvJac[0]= ( y_eta*z_zeta-y_zeta*z_eta ) *invdet;
    InvJac[3]=- ( x_eta*z_zeta-x_zeta*z_eta ) *invdet; //
    InvJac[6]= ( y_zeta*x_eta-x_zeta*y_eta ) *invdet;
    InvJac[1]=- ( y_xi*z_zeta-y_zeta*z_xi ) *invdet; //
    InvJac[4]= ( x_xi*z_zeta-x_zeta*z_xi ) *invdet;
    InvJac[7]=- ( y_zeta*x_xi-y_xi*x_zeta ) *invdet; //
    InvJac[2]= ( z_eta*y_xi-y_eta*z_xi ) *invdet;
    InvJac[5]=- ( x_xi*z_eta-x_eta*z_xi ) *invdet; //
    InvJac[8]= ( x_xi*y_eta-y_xi*x_eta ) *invdet;


    return ( det );
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac1D_nodes ( // 2D  jacobean ->
    const int ng,   // Gaussian point <-
    double x[],       // coordinates  <-
    double InvJac[] // Jacobean
)  {// ===================================================

    double x_xi=0.;
    int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss1[0];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[0]+ng;
        x_xi  +=x[s]*_dphidxez_map1_nodes[0][sng];
    }
    double det=x_xi;
    double idet=1./det;
    InvJac[0]=x_xi*idet;     //deta dy
    return ( det );
}

// =======================================================
/// This function computes the 2D derivatives and
/// Jacobian at the gaussian point ng
double MGFE::Jac2D_nodes ( // 2D  jacobean ->
    const int ng,   // Gaussian point <-
    double x[],       // coordinates  <-
    double InvJac[] // Jacobean
)  {// ===================================================

    double x_xi=0., x_eta =0.,y_xi =0., y_eta=0.;
    int nshape=_NoShape[1];
    int offset=_NoShape[1]*_NoGauss1[1];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[1]+ng;
        x_xi  +=x[s]*_dphidxez_map1_nodes[1][sng];
        x_eta +=x[s]*_dphidxez_map1_nodes[1][sng+offset];
        y_xi  +=x[s+NDOF_FEM]*_dphidxez_map1_nodes[1][sng];
        y_eta +=x[s+NDOF_FEM]*_dphidxez_map1_nodes[1][sng+offset];
    }

//   std::cout<<x_xi<<"  "<<x_eta<<"  "<<y_xi<<"  "<<y_eta<<std::endl;
//
    double det=x_xi*y_eta-y_xi*x_eta;
    double idet=1./det;
    InvJac[0]=y_eta*idet;    // dxi dx
    InvJac[1]=-y_xi*idet;    //deta dx
    InvJac[2]=-x_eta*idet;   // dxi dy
    InvJac[3]=x_xi*idet;     //deta dy

    return ( det );
}


// ==========================================
double MGFE::JacSur ( // surface jacobean ->
    const int ng,    // gaussian point
    double x[],       // coordinates
    double InvJac[]  // Jacobean ->
) const { // ================================

    double JacSur=0.;
#if  DIMENSION==1
    JacSur = JacSur1D ( ng,x,InvJac );
#endif
#if  DIMENSION==2
    JacSur = JacSur2D ( ng,x,InvJac );
#endif
#if  DIMENSION==3
    JacSur = JacSur3D ( ng,x,InvJac );
#endif
    return JacSur;
}
// =====================================
/// This function compute the line Jacobian at the gaussian point ng
double MGFE::JacSur1D ( // 1D surface jacobean ->
    const int ng,       // gaussian point  <-
    double x[],       // coordinates <-
    double InvJac[]  // Jacobean ->
) const { // ==========================

    // Values to compute at gaussian points
    double dxdxi=0.;//
//   double dydxi=0.;// d(x,y)d(xi,eta)
    const int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss[0];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[0]+ng;
        double dphidxi=_dphidxez_map1[0][sng];
        dxdxi += x[s]*dphidxi;
//     dydxi += x[s+NDOF_FEMB]*dphidxi;
    }
    //surface weighted jacobean
    double det=sqrt ( dxdxi*dxdxi );
    return ( det );
}

// =====================================
/// This function compute the line Jacobian at the gaussian point ng
double MGFE::JacSur2D ( // 2D surface jacobean ->
    const int ng,       // gaussian point  <-
    double x[],       // coordinates <-
    double InvJac[]  // Jacobean ->
) const { // ==========================

    // Values to compute at gaussian points
    double dxdxi=0.;//
    double dydxi=0.;// d(x,y)d(xi,eta)
    const int nshape=_NoShape[0];
//   int offset=_NoShape[0]*_NoGauss[0];
    for ( int s=0; s<nshape; s++ ) {
        int sng=s*_NoGauss1[0]+ng;
        double dphidxi=_dphidxez_map1[0][sng];
        dxdxi += x[s]*dphidxi;
        dydxi += x[s+NDOF_FEMB]*dphidxi;
    }
    //It functions only on straight boundary edge!!!!!!!!!!!!!!!!!!!!
    InvJac[0] = 1./ ( dxdxi+dydxi ); //porcata!!!!!!!!!!!
    //surface weighted jacobean
    double det=sqrt ( dxdxi*dxdxi+dydxi*dydxi );
    return ( det );
}



// ===============================================
/// This function computes the surface Jacobian at the gaussian point ng:
double MGFE::JacSur3D (
    const int ng,// gaussian point <-
    double x[],       // coordinates <-
    double InvJac[]  // Jacobean ->
) const {// ======================================

    // Values to compute at gaussian points
    double dxdxi=0.;
    double dxdeta=0.;//
    double dydxi=0.;
    double dydeta=0.;// d(x,y,z)d(xi,eta)
    double dzdxi=0.;
    double dzdeta=0.;//

//   int nshape=_NoShape[1];
    int offset=_NoShape[1]*_NoGauss1[1];

    for ( int s=0; s< ( int ) _NoShape[1]; s++ ) {

        int sng=s*_NoGauss1[1]+ng;
        double dphidxi=_dphidxez_map1[1][sng];
        double dphideta=_dphidxez_map1[1][sng+offset];

        dxdxi += x[s]*dphidxi;
        dxdeta += x[s]*dphideta;
        dydxi += x[s+NDOF_FEMB]*dphidxi;
        dydeta += x[s+NDOF_FEMB]*dphideta;
        dzdxi += x[s+2*NDOF_FEMB]*dphidxi;
        dzdeta += x[s+2*NDOF_FEMB]*dphideta;
    }
    //surface weighted jacobean
    double det=sqrt ( ( dxdxi*dydeta-dxdeta*dydxi ) * ( dxdxi*dydeta-dxdeta*dydxi ) +
    ( dydxi*dzdeta-dydeta*dzdxi ) * ( dydxi*dzdeta-dydeta*dzdxi ) +
    ( dzdxi*dxdeta-dzdeta*dxdxi ) * ( dzdxi*dxdeta-dzdeta*dxdxi )
                    );
    double idet = 1./det;
    InvJac[0]=dydeta*idet;    // dxi dx
    InvJac[1]=-dydxi*idet;    //deta dx
    InvJac[2]=-dxdeta*idet;   // dxi dy
    InvJac[3]=dxdxi*idet;     //deta dy
    return ( det );
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
void MGFE::get_phi_gl_g (
    const int kdim, // dimension <-
    const int qp,   // gaussian point <-
    double phi[]     // shape functions ->
) {// =================================================
    for ( int ish=0; ish<_NoShape[kdim-1]; ish++ ) {
        phi[ish] = _phi_map1[kdim-1][ish*_NoGauss1[kdim-1]+qp];
    }
    return;
}

// ==================================================================
/// This function computes the shape values at the gauss point qp
void MGFE::get_phi_gl_g (
    const int kdim,                       // dimension        <-
    const int qp,                         // gaussian point   <-
    std::vector<double> & phi             // shape functions  ->
) {// =====================================
    for ( int ish=0; ish<_NoShape[kdim-1]; ish++ ) {
        phi[ish] = _phi_map1[kdim-1][ish*_NoGauss1[kdim-1]+qp];
    }
    return;
}
void MGFE::get_phi_g_arb_el (
    const int kdim, // dimension <-
    const int qp,   // gaussian point <-
    double phi[],     // shape functions ->
    int FamilyType
) {// =================================================
    for ( int ish=0; ish<_GNoShape[FamilyType][kdim-1]; ish++ ) {
        phi[ish] = _Gphi_map1[ ( kdim-1 ) + FamilyType*3][ish*_GNoGauss1[FamilyType][kdim-1]+qp];
    }
    return;
}
// =================================================================
/// Shape functions derivatives dphi[id+i* el_nnodes] =Ti
/// Tx Ty Tz
///  tensor  order  at the gauss point qp, node id:
///  we have  dphi[Tx(0),Tx(1),..Tx(id), Ty(0),Ty(1),..Ty(id), Tz(0),Tz(1),..Tz(id) ]
void   MGFE::get_dphi_gl_g (
    const int kdim,// dimension <-
    const int qp,  // gaussian point <
    const double InvJac[], // Jacobean
    double dphi[]   // global derivatives ->
) {// =========================================

    const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
    const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
    const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset

    double gradphi_g[DIMENSION]; // temp grad phi

    for ( int eln=0; eln<el_nnodes; eln++ )    {
        int lqp=eln* el_ngauss+qp;
        for ( int idim=0; idim<_dim; idim++ ) gradphi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
        for ( int idim=0; idim<_dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<_dim; jdim++ ) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
            dphi[eln+idim* el_nnodes] = sum;
        }
    }
    return;
}


void   MGFE::get_dphi_g_arb_el (
    const int kdim,// dimension <-
    const int qp,  // gaussian point <
    const double InvJac[], // Jacobean
    double dphi[],   // global derivatives ->
    int FamilyType
) {// =========================================

    const int el_nnodes =  _GNoShape[FamilyType][kdim-1];             // # of shape functions
    const int el_ngauss =  _GNoGauss1[FamilyType][kdim-1];            // # of gauss points
    const int goffset   =  el_nnodes*_GNoGauss1[FamilyType][kdim-1];  //gauss offset

    double gradphi_g[DIMENSION]; // temp grad phi

    for ( int eln=0; eln<el_nnodes; eln++ )    {
        int lqp=eln* el_ngauss+qp;
        for ( int idim=0; idim<_dim; idim++ ) gradphi_g[idim] = _Gdphidxez_map1[kdim-1 + 3*FamilyType][lqp+idim*goffset];
        for ( int idim=0; idim<_dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<_dim; jdim++ ) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
            dphi[eln+idim* el_nnodes] = sum;
        }
    }


    return;
}
// =============================================================
/// Shape functions 2nd derivatives at the gauss point qp
/// ddphi[id*_dim*_dim+i*_dim+j]= Tij
///    Txx Txy Txz
///    Tyx Tyy Tyz
///    Tzx Tzy Tzz
///
///     ddphi[Txx(0),Txy(0),Txz(0),Tyx(0),Tyy(0),Tyz(0),Tzx(0),Tzy(0),Tzz(0),
///                 ....................................
///                 Txx(id),Txy(id),Txz,Tyx(id),Tyy,Tyz,Tzx(id),Tzy(id),Tzz(id)]
///   tensor  order  at the gauss point qp, node inode
void   MGFE::get_ddphi_gl_g (
    const int kdim,// dimension <-
    const int qp,  // gaussian point <
    const double InvJac[], // Jacobean
    double ddphi[]   // global derivatives ->
) {// =========================================

    const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
    const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
    const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset
    double ddphi_loc[NDOF_FEM*DIMENSION*DIMENSION];

    for ( int eln=0; eln<el_nnodes; eln++ ) { // LOOP OVER NODES ===================================
        const int shift = eln*_dim*_dim;   // offset for node derivatives
        const int lqp=eln*el_ngauss+qp;
        for ( int idim=0; idim<_dim*_dim; idim++ )
            ddphi_loc[shift+idim] = _dphidxx_map1[kdim-1][lqp+idim*goffset];


        double Hess_tmp[DIMENSION*DIMENSION],  JacTf[DIMENSION*DIMENSION];
        for ( int init=0; init<_dim; init++ ) for ( int jnit=0; jnit<_dim; jnit++ ) {
                Hess_tmp[init*_dim+jnit]=0.;
                ddphi[shift+init*_dim+jnit]=0.;
                JacTf[init*_dim+jnit] = InvJac[jnit*_dim+init];
            }

        for ( int ii=0; ii<_dim; ii++ ) for ( int jj=0; jj<_dim; jj++ )  for ( int ss=0; ss<_dim; ss++ )
                    Hess_tmp[ii*_dim+jj] += ddphi_loc[shift+ii*_dim+ss]*JacTf[ss*_dim+jj];

        // calculting d2 phi / dxi dxj
        for ( int ii=0; ii<_dim; ii++ )
            for ( int jj=0; jj<_dim; jj++ )
                for ( int ss=0; ss<_dim; ss++ )
                    ddphi[shift+ii*_dim+jj] += InvJac[ii*_dim+ss]*Hess_tmp[ss*_dim+jj];
    }// END LOOP OVER NODES =====================================================================

    //Pay attention, different order with respect to first order derivatives!
//   delete[]ddphi_loc;
    return;
}
// =============================================================
/// Shape functions derivatives at the nodal points
void   MGFE::get_dphi_node (
    const int kdim,// dimension <-
    const int node,  // nodal point <-
    const double InvJac[], // Jacobean <-
    double dphi[]   // global derivatives at the nodal point ->
) {// =========================================

    const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
    const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
    const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset

    double gradphi_g[DIMENSION]; // temp grad phi

    for ( int eln=0; eln<el_nnodes; eln++ )    {
        int lqp=eln* el_ngauss+node;
        for ( int idim=0; idim<_dim; idim++ ) gradphi_g[idim] = _dphidxez_map1_nodes[kdim-1][lqp+idim*goffset];
        for ( int idim=0; idim<_dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<_dim; jdim++ ) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
            dphi[eln+idim* el_nnodes] = sum;
        }
    }
    return;
}


// =============================================================
/// Shape functions derivatives at the nodal points
void   MGFE::get_dphi_arb_node ( std::vector<double> NodeCoord,
                                 const int order,
                                 double InvJac[],
double dphi[] ) { // =========================================
    const int kdim = NodeCoord.size();

    const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
    const int el_ngauss =  _NoGauss1[kdim-1];            // # of gauss points
    const int goffset   =  el_nnodes*_NoGauss1[kdim-1];  //gauss offset

    double gradphi_g[DIMENSION]; // temp grad phi

    int el_nodes = ( order==1 ) ? 4:9;
    if ( kdim==3 ) el_nodes *= ( order==1 ) ? 2:3;

    for ( int eln=0; eln<el_nodes; eln++ )    {
        gradphi_g[0] = 1.;
        gradphi_g[1] = 1.;
        gradphi_g[kdim-1] = 1.;
        for ( int idim=0; idim<kdim; idim++ ) { // derivative of eln-th phi ind idim direction
            const int idim2 = ( idim+1 ) %kdim;
            gradphi_g[idim] = 0.5*_CooQ9[eln+idim*_Q9Off] * ( 0.5* ( 1.+NodeCoord[idim2]*_CooQ9[eln+idim2*_Q9Off] ) );
        }
        for ( int idim=0; idim<_dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<_dim; jdim++ ) sum += InvJac[jdim+idim*_dim]*gradphi_g[jdim];
            dphi[eln+idim* el_nodes] = sum;
        }
    }
    return;
}

// =============================================================
/// Shape functions derivatives at the gauss point qp
void   MGFE::get_dphi_gl_g (
    const int kdim,// dimension <-
    const int qp,  // gaussian point <
    const double InvJac[], // Jacobean
    double dphi[],   // global derivatives ->
    int sdim
) {// =========================================

    const int el_nnodes =  _NoShape[kdim-1];             // # of shape functions
    const int el_ngauss =  _NoGauss1[kdim-1];            // # of guass points
    const int goffset   =  el_nnodes*el_ngauss;  //gauss offset

    double gradphi_g[DIMENSION]; // temp grad phi
    for ( int idim=0; idim<3; idim++ ) gradphi_g[idim] =0.;
    for ( int eln=0; eln<el_nnodes; eln++ )    {
        int lqp=eln* el_ngauss+qp;
        for ( int idim=0; idim<sdim; idim++ ) gradphi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
        for ( int idim=0; idim<sdim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<sdim; jdim++ ) sum += InvJac[jdim+idim*sdim]*gradphi_g[jdim];
            dphi[eln+idim* el_nnodes] = sum;
        }
    }
    return;
}

// =============================================================
/// Shape functions derivatives at the gauss point qp
void   MGFE::get_dphi_gl_g (
    const int kdim,           // dimension <-
    const int qp,             // gaussian point <-
    const double InvJac[], // Jacobean
    std::vector<double> & dphi // global derivatives ->
) {// =========================================

    const int el_nnodes =    _NoShape[kdim-1];           // # of shape functions
    const int el_ngauss   =  _NoGauss1[kdim-1];          // # of guass points
    const int goffset   =    el_nnodes*_NoGauss1[kdim-1];   //gauss offset

    double dphidxi_g[DIMENSION];// temp grad phi

    for ( int eln=0; eln<el_nnodes; eln++ )    {
        int lqp=eln* el_ngauss+qp;
        for ( int idim=0; idim<_dim; idim++ ) dphidxi_g[idim] = _dphidxez_map1[kdim-1][lqp+idim*goffset];
        for ( int idim=0; idim<_dim; idim++ ) {
            double sum = 0.;
            for ( int jdim=0; jdim<_dim; jdim++ ) sum += InvJac[jdim+idim*_dim]*dphidxi_g[jdim];
            dphi[eln+idim* el_nnodes] = sum;
        }
    }
    return;
}

// ======================================
///This function computes the normal at the gauss point
void MGFE::normal_g (
    const double* xx,   // all surface coordinates <-
    double* normal_g    // normal ->
) const {// ======================================

// coordinates
    double xx3D[3*NDOF_FEMB];
    for ( int i=0; i<_dim*NDOF_FEMB; i++ ) xx3D[i]=xx[i];

    double tg01[DIMENSION];  // tangent  line
#if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
    for ( int i=0; i<2; i++ )  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
    normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
    normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
#else // 3D ---------------------------------------------------
    double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
    for ( int i=0; i<3; i++ )  {
        tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
        tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
        tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
    }
//   _mgutils.cross(tg01,tg03,normal_g);
    normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
    normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
    normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];

#endif
// normalization -----------------------
    double mm=0.;
    for ( int idim=0; idim< _dim; idim++ ) mm +=normal_g[idim]*normal_g[idim];
    mm=sqrt ( mm );
    for ( int idim=0; idim< _dim; idim++ ) normal_g[idim] /=mm;

    return;

}


// ======================================
///This function computes the normal at the gauss point
///   x_c[] is an interior point
void MGFE::normal_g (
    const double* xx,   // all surface coordinates <-
    const double x_c[],   // point inside the element <-
    double* normal_g    // normal ->
) const {// ======================================

// coordinates
    double xx3D[3*NDOF_FEMB];
    for ( int i=0; i<_dim*NDOF_FEMB; i++ ) xx3D[i]=xx[i];

    double tg01[DIMENSION];  // tangent  line
#if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
    for ( int i=0; i<2; i++ )  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
    normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
    normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
    // if the sign is not correct then reverse it

    if ( normal_g[0]* ( x_c[0]-xx3D[0] ) +normal_g[1]* ( x_c[1]-xx3D[0+NDOF_FEMB] ) >0. ) {
        normal_g[0] *=-1.;
        normal_g[1] *=-1.;
//    std::cout << " Normal inverted ! ------------------------------------  \n";
    }

#else // 3D ---------------------------------------------------
    double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
    for ( int i=0; i<3; i++ )  {
        tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
        tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
        tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
    }
//   _mgutils.cross(tg01,tg03,normal_g);
    normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
    normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
    normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];

    // if the sign is not correct then reverse it
    if ( normal_g[0]* ( x_c[0]-xx3D[0] ) +
    normal_g[1]* ( x_c[1]-xx3D[0+NDOF_FEMB] ) +
    normal_g[2]* ( x_c[2]-xx3D[0+2*NDOF_FEMB] ) >0 ) {
        normal_g[0] *=-1;
        normal_g[1] *=-1;
        normal_g[2] *=-1;
    }

#endif
// normalization -----------------------
    double mm=0.;
    for ( int idim=0; idim< _dim; idim++ ) mm +=normal_g[idim]*normal_g[idim];
    mm=sqrt ( mm );
    for ( int idim=0; idim< _dim; idim++ ) normal_g[idim] /=mm;

    return;

}


// ======================================
///This function computes the normal at the gauss point
///   x_c[] is an interior point
void MGFE::normal_g (
    const double* xx,   // all surface coordinates <-
    const double x_c[],   // point inside the element <-
    double* normal_g,    // normal ->
    int & sign_normal
) const {// ======================================

// coordinates
    double xx3D[3*NDOF_FEMB];
    for ( int i=0; i<_dim*NDOF_FEMB; i++ ) xx3D[i]=xx[i];

    double tg01[DIMENSION];  // tangent  line
#if  DIMENSION==2 //  2D -----------------------------------------------
//  The surface elements are such that when you go from the 1st to
//   the 2nd point, the outward normal is to the RIGHT
    for ( int i=0; i<2; i++ )  tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
//  Rotation matrix (0 1; -1 0), 90deg clockwise
    normal_g[0] =  0.*tg01[0] + 1.*tg01[1];
    normal_g[1] = -1.*tg01[0] + 0.*tg01[1];
    // if the sign is not correct then reverse it
    sign_normal=1;
    if ( normal_g[0]* ( x_c[0]-xx3D[0] ) +normal_g[1]* ( x_c[1]-xx3D[0+NDOF_FEMB] ) >0. ) {
        normal_g[0] *=-1.;
        normal_g[1] *=-1.;
        sign_normal=-1;
//    std::cout << " Normal inverted ! ------------------------------------  \n";
    }

#else // 3D ---------------------------------------------------
    double tg03[DIMENSION];  // tangent plane
//   the cross product of the two tangent vectors
    for ( int i=0; i<3; i++ )  {
        tg01[i]= xx3D[1+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#if ELTYPE == 27
        tg03[i]= xx3D[3+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
#if ELTYPE == 10 //depend of the orentation of points
        tg03[i]= xx3D[2+i*NDOF_FEMB] - xx3D[0+i*NDOF_FEMB];
#endif
    }
//   _mgutils.cross(tg01,tg03,normal_g);
    normal_g[0]=tg01[1]*tg03[2]-tg01[2]*tg03[1];
    normal_g[1]=tg01[2]*tg03[0]-tg01[0]*tg03[2];
    normal_g[2]=tg01[0]*tg03[1]-tg01[1]*tg03[0];

    // if the sign is not correct then reverse it
    if ( normal_g[0]* ( x_c[0]-xx3D[0] ) +
    normal_g[1]* ( x_c[1]-xx3D[0+NDOF_FEMB] ) +
    normal_g[2]* ( x_c[2]-xx3D[0+2*NDOF_FEMB] ) >0 ) {
        normal_g[0] *=-1;
        normal_g[1] *=-1;
        normal_g[2] *=-1;
    }

#endif
// normalization -----------------------
    double mm=0.;
    for ( int idim=0; idim< _dim; idim++ ) mm +=normal_g[idim]*normal_g[idim];
    mm=sqrt ( mm );
    for ( int idim=0; idim< _dim; idim++ ) normal_g[idim] /=mm;

    return;

}


int MGFE::GetFamilyType ( int elem_dof, int dim ) {
    int FamilyType;
    switch ( dim ) {
    case 2:
        if ( elem_dof == 3 || elem_dof == 6 ) FamilyType = 0;
        else if ( elem_dof == 4 || elem_dof == 9 ) FamilyType = 1;
        else printf ( "MGFE::GetFamilyType Unknown FamilyType for element with dimension %d and %d number of dof \n",dim,elem_dof );
        break;
    case 3:
        if ( elem_dof == 7 || elem_dof == 10 ) FamilyType = 0;
        else if ( elem_dof == 8 || elem_dof == 27 ) FamilyType = 1;
        else printf ( "MGFE::GetFamilyType Unknown FamilyType for element with dimension %d and %d number of dof \n",dim,elem_dof );
        break;
    default :
        FamilyType = 1;
        break;
    }
    return FamilyType;
}

double MGFE::Tri_2d_LinearPhi ( int nPhi, double point[] ) {
    double PhiVal;
    int a = _CooTriEl[nPhi*3];
    int b = _CooTriEl[nPhi*3 +1];
    int c = _CooTriEl[nPhi*3 +2];
    PhiVal = a*point[0] + b*point[1] + c;
    return PhiVal;
}

 
double MGFE::Tri_2d_QuadraticPhi ( int nPhi, double point[] ) {
    double PhiVal;
    int m = nPhi/3;
    int l1 = ( nPhi%3 );
    int l2 = ( l1+1 ) %3;
    double L1, L2;
    L1 = Tri_2d_LinearPhi ( l1,point );
    L2 = Tri_2d_LinearPhi ( l2,point );
    PhiVal = L1* ( ( 1-m ) * ( 2.*L1 -1 ) + m*4.*L2 );
    return PhiVal;
}

double MGFE::Tri_2d_LinearDerPhi ( int nPhi, double point[], int dir ) {
    double PhiDer;
    PhiDer = _CooTriEl[nPhi*3 + dir];
    return PhiDer;
}

double MGFE::Tri_2d_QuadraticDerPhi ( int nPhi, double point[], int dir ) {
    double PhiDer;
    int m = nPhi/3;
    int l1 = ( nPhi%3 );
    int l2 = ( l1+1 ) %3;
    double L1, L2, dL1, dL2;
    L1 = Tri_2d_LinearPhi ( l1,point );
    L2 = Tri_2d_LinearPhi ( l2,point );
    dL1 = Tri_2d_LinearDerPhi ( l1, point, dir );
    dL2 = Tri_2d_LinearDerPhi ( l2, point, dir );

    PhiDer = dL1* ( ( 1-m ) * ( 2.*L1 -1 ) + m*4.*L2 ) + L1* ( ( 1-m ) * ( 2.*dL1 ) + m*4.*dL2 );
    return PhiDer;
}

double MGFE::Tri_2d_QuadraticDer2Phi ( int nPhi, double point[], int dir1, int dir2 ) {
    double Phi2Der;
    int m = nPhi/3;
    int l1 = ( nPhi%3 );
    int l2 = ( l1+1 ) %3;
    double dL1_2, dL2_2, dL1_1, dL2_1;

    dL1_1 = Tri_2d_LinearDerPhi ( l1, point, dir1 );
    dL2_1 = Tri_2d_LinearDerPhi ( l2, point, dir1 );
    dL1_2 = Tri_2d_LinearDerPhi ( l1, point, dir2 );
    dL2_2 = Tri_2d_LinearDerPhi ( l2, point, dir2 );

    Phi2Der = dL1_1* ( ( 1-m ) * ( 2.*dL1_2 ) + m*4.*dL2_2 ) + dL1_2* ( ( 1-m ) * ( 2.*dL1_1 ) + m*4.*dL2_1 );
    return Phi2Der;
}



double MGFE::Edge_Quad_Phi ( int PhiCoeff, double Coordinate ) {
    double PhiVal;
    PhiVal = ( 1-0.5*fabs ( PhiCoeff ) ) * ( ( 2*fabs ( PhiCoeff )-1 ) *Coordinate*Coordinate + PhiCoeff*Coordinate + ( 1-fabs ( PhiCoeff ) ) );
    return PhiVal;
}

double MGFE::Edge_Quad_DPhi ( int PhiCoeff, double Coordinate ) {
    double DPhiVal;
    DPhiVal = ( 1-0.5*fabs ( PhiCoeff ) ) * ( 2.* ( 2*fabs ( PhiCoeff )-1 ) *Coordinate + PhiCoeff );
    return DPhiVal;
}

double MGFE::Edge_Quad_D2Phi ( int PhiCoeff, double Coordinate ) {
    double D2PhiVal;
    D2PhiVal = ( 1-0.5*fabs ( PhiCoeff ) ) * ( 2.* ( 2*fabs ( PhiCoeff )-1 ) );
    return D2PhiVal;
}

double MGFE::Edge_Lin_Phi ( int PhiCoeff, double Coordinate ) {
    double PhiVal;
    PhiVal = 0.5* ( 1+ PhiCoeff*Coordinate );
    return PhiVal;
}

double MGFE::Edge_Lin_DPhi ( int PhiCoeff, double Coordinate ) {
    double DPhiVal;
    DPhiVal = 0.5*PhiCoeff;
    return DPhiVal;
}

double MGFE::Rec_Lin_Phi ( int nPhi, double point[], int dimension ) {
    double PhiVal = 1.;
    for ( int dir=0; dir<dimension; dir++ ) {
        int lambda_i;
        if ( dimension==1 ) lambda_i = _CooE3[nPhi];
        else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + dir*_Q9Off];
        else if ( dimension==3 ) lambda_i = _CooH27[nPhi + dir*_H27Off];
        else printf ( "\033[1;31m MGFE::Rec_Lin_Phi unkown phi for dimension %d \n \033[0m",dimension );

        PhiVal *= Edge_Lin_Phi ( lambda_i, point[dir] );
    }
    return PhiVal;
}

double MGFE::Rec_Quad_Phi ( int nPhi, double point[], int dimension ) {
    double PhiVal = 1.;
    for ( int dir=0; dir<dimension; dir++ ) {
        int lambda_i;
        if ( dimension==1 ) lambda_i = _CooE3[nPhi];
        else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + dir*_Q9Off];
        else if ( dimension==3 ) lambda_i = _CooH27[nPhi + dir*_H27Off];
        else printf ( "\033[1;31m MGFE::Rec_Quad_Phi unkown phi for dimension %d \n \033[0m",dimension );

        PhiVal *= Edge_Quad_Phi ( lambda_i, point[dir] );
    }
    return PhiVal;
}

double MGFE::Rec_Lin_DPhi ( int nPhi, double point[], int dimension, int DirDer ) {
    int lambda_i;
    if ( dimension==1 ) lambda_i = _CooE3[nPhi];
    else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + DirDer*_Q9Off];
    else if ( dimension==3 ) lambda_i = _CooH27[nPhi + DirDer*_H27Off];
    else printf ( "\033[1;31m MGFE::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m",dimension );

    double PhiDer = Edge_Lin_DPhi ( lambda_i, point[DirDer] );

    for ( int dir=DirDer+1; dir<dimension+DirDer; dir++ ) {
        int lambda_i;
        int direction = dir%dimension;
        if ( dimension==1 ) lambda_i = _CooE3[nPhi];
        else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + direction*_Q9Off];
        else if ( dimension==3 ) lambda_i = _CooH27[nPhi + direction*_H27Off];
        else printf ( "\033[1;31m MGFE::Rec_Lin_DPhi unkown phi for dimension %d \n \033[0m",dimension );
        PhiDer *= Edge_Lin_Phi ( lambda_i, point[direction] );
    }
    return PhiDer;
}

double MGFE::Rec_Quad_DPhi ( int nPhi, double point[], int dimension, int DirDer ) {
    int lambda_i;
    if ( dimension==1 ) lambda_i = _CooE3[nPhi];
    else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + DirDer*_Q9Off];
    else if ( dimension==3 ) lambda_i = _CooH27[nPhi + DirDer*_H27Off];
    else printf ( "\033[1;31m MGFE::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m",dimension );

    double PhiDer = Edge_Quad_DPhi ( lambda_i, point[DirDer] );

    for ( int dir=DirDer+1; dir<dimension+DirDer; dir++ ) {
        int lambda_i;
        int direction = dir%dimension;
        if ( dimension==1 ) lambda_i = _CooE3[nPhi];
        else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + direction*_Q9Off];
        else if ( dimension==3 ) lambda_i = _CooH27[nPhi + direction*_H27Off];
        else printf ( "\033[1;31m MGFE::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m",dimension );
        PhiDer *= Edge_Quad_Phi ( lambda_i, point[direction] );
    }
    return PhiDer;
}

double MGFE::Rec_Lin_D2Phi ( int nPhi, double point[], int dimension, int DirDer1, int DirDer2 ) {
    double PhiDer;
    if ( DirDer1==DirDer2 ) PhiDer = 0.;
    else {
        PhiDer = 1;
        for ( int dir=0; dir<dimension; dir++ ) {
            int lambda_i;
            int direction = dir%dimension;
            if ( dimension==1 ) lambda_i = _CooE3[nPhi];
            else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + direction*_Q9Off];
            else if ( dimension==3 ) lambda_i = _CooH27[nPhi + direction*_H27Off];
            else printf ( "\033[1;31m MGFE::Rec_Lin_D2Phi unkown phi for dimension %d \n \033[0m",dimension );
            if ( dir==DirDer1 || dir==DirDer2 ) PhiDer *= Edge_Lin_DPhi ( lambda_i, point[dir] );
            else PhiDer *= Edge_Lin_Phi ( lambda_i, point[dir] );
        }
    }
    return PhiDer;
}

double MGFE::Rec_Quad_D2Phi ( int nPhi, double point[], int dimension, int DirDer1, int DirDer2 ) {
    double PhiDer;
    if ( DirDer1==DirDer2 ) {
        int lambda_i;
        if ( dimension==1 ) lambda_i = _CooE3[nPhi];
        else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + DirDer1*_Q9Off];
        else if ( dimension==3 ) lambda_i = _CooH27[nPhi + DirDer1*_H27Off];
        PhiDer = Edge_Quad_D2Phi ( lambda_i, point[DirDer1] );
        for ( int dir=DirDer1+1; dir<dimension+DirDer1; dir++ ) {
            int lambda_i;
            int direction = dir%dimension;
            if ( dimension==1 ) lambda_i = _CooE3[nPhi];
            else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + direction*_Q9Off];
            else if ( dimension==3 ) lambda_i = _CooH27[nPhi + direction*_H27Off];
            else printf ( "\033[1;31m MGFE::Rec_Quad_DPhi unkown phi for dimension %d \n \033[0m",dimension );
            PhiDer *= Edge_Quad_Phi ( lambda_i, point[direction] );
        }
    } else {
        PhiDer = 1;
        for ( int dir=0; dir<dimension; dir++ ) {
            int lambda_i;
            if ( dimension==1 ) lambda_i = _CooE3[nPhi];
            else if ( dimension==2 ) lambda_i = _CooQ9[nPhi + dir*_Q9Off];
            else if ( dimension==3 ) lambda_i = _CooH27[nPhi + dir*_H27Off];
            else printf ( "\033[1;31m MGFE::Rec_Quad_D2Phi unkown phi for dimension %d \n \033[0m",dimension );
            if ( dir==DirDer1 || dir==DirDer2 ) PhiDer *= Edge_Quad_DPhi ( lambda_i, point[dir] );
            else PhiDer *= Edge_Quad_Phi ( lambda_i, point[dir] );
        }
    }
    return PhiDer;
}

double MGFE::Tri_3d_LinearPhi ( int nPhi, double point[] ) {
    double PhiVal;
    int a = _CooTetra4[nPhi*4];
    int b = _CooTetra4[nPhi*4 +1];
    int c = _CooTetra4[nPhi*4 +2];
    int d = _CooTetra4[nPhi*4 +3];
    PhiVal = a*point[0] + b*point[1] + c*point[2] + d;
    return PhiVal;
}

double MGFE::Tri_3d_LinearDerPhi ( int nPhi, double point[], int dir ) {
    double PhiDer;
    PhiDer = _CooTetra4[nPhi*4 + dir];
    return PhiDer;
}

double MGFE::Tri_3d_QuadraticPhi ( int nPhi, double point[] ) {
    double PhiVal;
    int l1 = _CooTetra10[nPhi*4];
    int l2 = _CooTetra10[nPhi*4 +1];
    int a  = _CooTetra10[nPhi*4 +2];
    int b  = _CooTetra10[nPhi*4 +3];
    double L1, L2;
    L1 = Tri_3d_LinearPhi ( l1,point );
    L2 = Tri_3d_LinearPhi ( l2,point );
    PhiVal = L1* ( a*L2 + b );
    return PhiVal;
}

double MGFE::Tri_3d_QuadraticDerPhi ( int nPhi, double point[], int dir ) {
    double PhiDer;
    int l1 = _CooTetra10[nPhi*4];
    int l2 = _CooTetra10[nPhi*4 +1];
    int a  = _CooTetra10[nPhi*4 +2];
    int b  = _CooTetra10[nPhi*4 +3];
    double L1, L2, dL1, dL2;
    L1 = Tri_3d_LinearPhi ( l1,point );
    L2 = Tri_3d_LinearPhi ( l2,point );
    dL1 = Tri_3d_LinearDerPhi ( l1,point, dir );
    dL2 = Tri_3d_LinearDerPhi ( l2,point, dir );
    PhiDer = dL1* ( a*L2 + b ) + L1*a*dL2;
    return PhiDer;
}

double MGFE::Tri_3d_QuadraticDer2Phi ( int nPhi, double point[], int dir1, int dir2 ) {
    double PhiDer;
    int l1 = _CooTetra10[nPhi*4];
    int l2 = _CooTetra10[nPhi*4 +1];
    int a  = _CooTetra10[nPhi*4 +2];
    double dL1_1, dL2_1, dL1_2, dL2_2;
    dL1_1 = Tri_3d_LinearDerPhi ( l1,point, dir1 );
    dL2_1 = Tri_3d_LinearDerPhi ( l2,point, dir1 );
    dL1_2 = Tri_3d_LinearDerPhi ( l1,point, dir2 );
    dL2_2 = Tri_3d_LinearDerPhi ( l2,point, dir2 );
    PhiDer = a* ( dL1_1*dL2_2 + dL1_2*dL2_1 );
    return PhiDer;
}

double MGFE::FirstDerivateOfLocalPhi ( int nPhi, double point[], int Dimension, int DirDer, int FamilyType ) {
    double PhiDer;
    if ( FamilyType==0 ) { // TRIANGULAR ELEMENTS
        if ( _order==1 ) {
            if ( Dimension==1 ) PhiDer = Rec_Lin_DPhi ( nPhi, point, Dimension, DirDer );
            else if ( Dimension==2 ) PhiDer = Tri_2d_LinearDerPhi ( nPhi, point, DirDer );
            else if ( Dimension==3 ) PhiDer = Tri_3d_LinearDerPhi ( nPhi, point, DirDer );
        }
        if ( _order==2 ) {
            if ( Dimension==1 ) PhiDer = Rec_Quad_DPhi ( nPhi, point, Dimension, DirDer );
            else if ( Dimension==2 ) PhiDer = Tri_2d_QuadraticDerPhi ( nPhi, point, DirDer );
            else if ( Dimension==3 ) PhiDer = Tri_3d_QuadraticDerPhi ( nPhi, point, DirDer );
        }
    }
    if ( FamilyType==1 ) { // QUADRANGULAR ELEMENTS
        if ( _order==1 ) PhiDer = Rec_Lin_DPhi ( nPhi, point, Dimension, DirDer );
        if ( _order==2 ) PhiDer = Rec_Quad_DPhi ( nPhi, point, Dimension, DirDer );
    }
    return PhiDer;
}

double MGFE::ComputeInverseMatrix ( double Matrix[], double InvMatrix[], int dim ) {
    double det=0.;
    if ( dim==3 ) {
        for ( int dir=0; dir<dim; dir++ ) {
            int sign = 1 - 2* ( dir%2 );
            int idx1 = ( dir +1 ) %dim;
            int idx2 = ( dir +2 ) %dim;
            det += sign* ( Matrix[idx1*dim + idx1]*Matrix[idx2*dim + idx2] - Matrix[idx1*dim + idx2]*Matrix[idx2*dim + idx1] );
        }
    }
    else if ( dim==2 ) {
        det = ( Matrix[0]*Matrix[3] - Matrix[1]*Matrix[2] );
    }
    else if ( dim==1 ) det = Matrix[0];

    double idet=1./det;

    if ( dim==1 ) InvMatrix[0]=Matrix[0]*idet;     //
    if ( dim==2 ) {
        InvMatrix[0] =  Matrix[3]*idet;    // dxi dx
        InvMatrix[1] = -Matrix[2]*idet;    // deta dx
        InvMatrix[2] = -Matrix[1]*idet;    // dxi dy
        InvMatrix[3] =  Matrix[0]*idet;    // deta dy
    }
    if ( dim==3 ) {
        for ( int row=0; row<dim; row++ ) {// LOOP OVER ROWS - row1<row2
            int sign_row =  1 - 2* ( row%2 );
            int row1 = ( ( row +1 ) %dim < ( row +2 ) %dim ) ? ( row +1 ) %dim : ( row +2 ) %dim;
            int row2 = ( ( row +1 ) %dim > ( row +2 ) %dim ) ? ( row +1 ) %dim : ( row +2 ) %dim;
            for ( int col=0; col<dim; col++ ) {// LOOP OVER COLUMNS  - col1<col2
                int col1 = ( ( col +1 ) %dim < ( col +2 ) %dim ) ? ( col +1 ) %dim : ( col +2 ) %dim;
                int col2 = ( ( col +1 ) %dim > ( col +2 ) %dim ) ? ( col +1 ) %dim : ( col +2 ) %dim;
                int sign_col =  1 - 2* ( col%2 );
                InvMatrix[row*dim + col] = sign_col*sign_row * idet *
                                           ( Matrix[row1*dim + col1]*Matrix[row2*dim + col2] - Matrix[row1*dim + col2]*Matrix[row2*dim + col1] );
            }
        }
    }
    return det;
}
