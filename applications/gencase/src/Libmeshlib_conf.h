#ifndef __libmeshlibconf_h__
#define __libmeshlibconf_h__







// ****************************************
//  LIBMESH
// *****************************************
// *************************************************
//This include file rules the Libmesh dependencies
// ************************************************


// #define LM_INIT   //TODO useless now  *********** LIB-TYPE dependency

// #define LM_REFCOUNT  //uncoupled         *********** LIB-TYPE dependency

//#define LM_REAL      //uncoupled          *********** INCLUDE-TYPE dependency


//***********who cares************
#define LM_GENCASE   //hopeless, for now  *********** LIB-TYPE dependency
// #define LM_PERFLOG  //forgotten        *********** LIB-TYPE dependency
//***********end who cares********



//*****************
#ifdef LM_REFCOUNT

//refcount needs init (and also debug mode, by the way)
    #ifndef LM_INIT   
      #define LM_INIT
    #endif

//also,if refcount starts then LM_REAL must also start, otherwise
                      //you get ambiguous references to Real
    #ifndef LM_REAL 
    #define LM_REAL
    #endif
    
  #include "reference_counted_object.h"

#endif
 //******************

//you may want to use LM_INIT either directly or because you use LM_REFCOUNT










#endif