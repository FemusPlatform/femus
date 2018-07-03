
 #ifndef equation_fields212
 #define equation_fields212

#include <map>
#include <vector>
 
 enum  FIELDS{
    MG_NavierStokes = 0,
    NS_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    NSZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    MG_FluidStructure =0,
    FS_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    FSZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    MG_StructuralMechanics=0,
    SM_F  =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMX_F =0,     // [0] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMY_F =1,     // [1] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    SMZ_F =2,     // [2] -> Navier-Stokes or FSI or SM (quadratic (2),NS_EQUATIONS)
    MG_Pressure=3,
    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
    MG_Temperature=4,
    T_F   =4,     // [4] -> Temperature   (quadratic (2),T_EQUATIONS)
    MG_DynamicalTurbulence=5,
    K_F   =5,     // [5] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EW_F  =6,     // [6] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    MG_ThermalTurbulence =7,
    KTT_F =7,     // [7] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EWTT_F =8,    // [8] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    MG_Displacement=9,
    SDS_F =9,    // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
    SDSX_F =9,    // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
    SDSY_F =10,   // [10]-> Displacement (quadratic (2), DS_EQUATIONS)
    SDSZ_F =11,   // [11]-> Displacement (quadratic (2), DS_EQUATIONS)
    MG_DA=12,
    DA_F   =12,   // [12]-> DA solver (quadratic (2), DA_EQUATIONS)
    // adjoint
    MG_AdjointNavierStoke=13,
    NSA_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAX_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAY_F  =14,   // [15]-> adjoint NS or FSI equations (Dimension)
    NSAZ_F  =15,   // [15]-> adjoint NS or FSI equations (Dimension)
    MG_AdjointFluidStructure=13,
    FSA_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAX_F  =13,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAY_F  =14,   // [15]-> adjoint NS or FSI equations (Dimension)
    FSAZ_F  =15,   // [15]-> adjoint NS or FSI equations (Dimension)
    MG_AdjointDA=16,
    DA_P   =16,   // [13]-> DA solver (piecewise, DA_EQUATIONS)
    MG_AdjointTemperature=17,
    TA_F   =17,   // [14]-> Temp adjoint
    MG_AdjointTurbulence =18,
    KA_F   =18,   // [18] + 19 Adjoint turbulence
    
    //control
    MG_ColorFunction=20,
    MG_Laplacian=20,
    CO_F   =20,    // [16]-> Color function for FSI equations
    CTRL_F  =21,    // [16]-> Color function for FSI equations
    CTRLX_F  =21,    // [16]-> Color function for FSI equations
    CTRLY_F  =22,    // [16]-> Color function for FSI equations
    CTRLZ_F  =23,    // [16]-> Color function for FSI equations
  };
  
  class FIELDS_class{
  public: 
       std::vector<FIELDS> myproblemP;
       std::map<std::string,FIELDS> map_str2field;
     
       FIELDS_class(){
    // equation
           
    map_str2field["MG_NavierStokes"] = NS_F;
    map_str2field["NS_F"]=NS_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["NSX_F"]=NSX_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["NSY_F"]=NSY_F;    // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["NSZ_F"]=NSZ_F;     // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["MG_FluidStructure"] =FS_F;
    map_str2field["FS_F"] =FS_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["FSX_F"]=FSX_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["FSY_F"]=FSY_F;     // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["FSZ_F"] =FSZ_F;     // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["MG_StructuralMechanics"]=SM_F;
    map_str2field["SM_F"]  =SM_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["SMX_F"] =SMX_F;     // [0] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["SMY_F"] =SMY_F;     // [1] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["SMZ_F"] =SMZ_F;     // [2] -> Navier-Stokes or FSI or SM (quadratic (2);NS_EQUATIONS)
    map_str2field["MG_Pressure"]=P_F;
    map_str2field["P_F"]   =P_F;     // [3] -> Pressure (linear (1);NS_EQUATIONS==0 or 2)
    map_str2field["MG_Temperature"]=T_F;
    map_str2field["T_F"]   =T_F;     // [4] -> Temperature   (quadratic (2);T_EQUATIONS)
    map_str2field["MG_DynamicalTurbulence"]=K_F;
    map_str2field["K_F"]   =K_F;     // [5] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
    map_str2field["EW_F"]  =EW_F;     // [6] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
    map_str2field["MG_ThermalTurbulence"] =KTT_F;
    map_str2field["KTT_F"] =KTT_F;     // [7] -> Turbulence K  (quadratic (2);TB K_EQUATIONS)
    map_str2field["EWTT_F"] =EWTT_F;    // [8] -> Turbulence W  (quadratic (2);TB W_EQUATIONS)
    map_str2field["MG_Displacement"]=SDS_F;
    map_str2field["SDS_F"] =SDS_F;    // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
    map_str2field["SDSX_F"] =SDSX_F;    // [9] -> Displacement (quadratic (2); DS_EQUATIONS)
    map_str2field["SDSY_F"] =SDSY_F;   // [10]-> Displacement (quadratic (2); DS_EQUATIONS)
    map_str2field["SDSZ_F"] =SDSZ_F;   // [11]-> Displacement (quadratic (2); DS_EQUATIONS)
    map_str2field["MG_DA"]=DA_F;
    map_str2field["DA_F"]   =DA_F;   // [12]-> DA solver (quadratic (2); DA_EQUATIONS)
    // adjoint
    map_str2field["MG_AdjointNavierStokes"]=NSA_F;
    map_str2field["NSA_F"]  =NSA_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["NSAX_F"]  =NSAX_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["NSAY_F"]  =NSAY_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["NSAZ_F"]  =NSAZ_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["MG_AdjointFluidStructure"]=FSA_F;
    map_str2field["FSA_F"]  =FSA_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["FSAX_F"]  =FSAX_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["FSAY_F"]  =FSAY_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["FSAZ_F"]  =FSAZ_F;   // [15]-> adjoint NS or FSI equations (Dimension)
    map_str2field["MG_AdjointDA"]=DA_P;
    map_str2field["DA_P"]   =DA_P;   // [13]-> DA solver (piecewise; DA_EQUATIONS)
    map_str2field["MG_AdjointTemperature"]=TA_F;
    map_str2field["TA_F"]   =TA_F;   // [14]-> Temp adjoint
    map_str2field["MG_AdjointTurbulence"] =KA_F;
    map_str2field["KA_F"]   =KA_F;   // [18] + 19 Adjoint turbulence

    //control
    map_str2field["MG_ColorFunction"]=CO_F;
    map_str2field["MG_Laplacian"]=CO_F;
    map_str2field["CO_F"]   =CO_F;    // [16]-> Color function for FSI equations
    map_str2field["CTRL_F"]  =CTRL_F;    // [16]-> Color function for FSI equations
    map_str2field["CTRLX_F"]  =CTRLX_F;    // [16]-> Color function for FSI equations
    map_str2field["CTRLY_F"]  =CTRLY_F;    // [16]-> Color function for FSI equations
    map_str2field["CTRLZ_F"]  =CTRLZ_F;    // [16]-> Color function for FSI equations
}
      
      
  void init_FIELDS(int iname);
 
      
};   

 #endif  // end file _equationsconf_
