 #ifndef _equation_fields_
 #define _equation_fields_

 
 enum  FIELDS{
    NS_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FS_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SM_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
    T_F   =4,     // [4] -> Temperature   (quadratic (2),T_EQUATIONS)
    K_F   =5,     // [5] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EW_F  =6,     // [6] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    KTT_F =7,     // [7] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EWTT_F =8,    // [8] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    SDS_F =9,    // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
    SDSX_F =9,    // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
    SDSY_F =10,   // [10]-> Displacement (quadratic (2), DS_EQUATIONS)
    SDSZ_F =11,   // [11]-> Displacement (quadratic (2), DS_EQUATIONS)
    DA_F   =12,   // [12]-> DA solver (quadratic (2), DA_EQUATIONS)
    // adjoint
    NSA_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAX_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAY_F  =14,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAZ_F  =15,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSA_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAX_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAY_F  =14,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAZ_F  =15,   // [15]-> adjoint NS or FSI equations (Dimension)
    DA_P   =16,   // [13]-> DA solver (piecewise, DA_EQUATIONS)
    TA_F   =17,   // [14]-> Temp adjoint
    KA_F   =18,   // [18] + 19 Adjoint turbulence
    
    //control
    CO_F   =20,    // [16]-> Color function for FSI equations
    CTRL_F  =21,    // [16]-> Color function for FSI equations
    CTRLX_F  =21,    // [16]-> Color function for FSI equations
    CTRLY_F  =22,    // [16]-> Color function for FSI equations
    CTRLZ_F  =23,    // [16]-> Color function for FSI equations
  };

 #endif  // end file _equationsconf_
