#ifndef _equationsconfw
#define _equationsconfw

 #define FSI_EQUATIONS (1)

 #ifdef FSI_EQUATIONS
 /// Displacement
 #define  DS_EQUATIONS
 #if (FSI_EQUATIONS%2==0)       // projection method
  #define FSIP_EQUATIONS (1)     // need separated P equations
 #endif
 #endif


#endif  // end file _equationsconf_
