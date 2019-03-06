#ifndef _printinfo_cfg
#define _printinfo_cfg

// ===========================
#define PRINT_INFO (1)
// PRINT_INFO (option=0) routine info

#define PRINT_TIME (1)
// PRINT_INFO (option=1) time info
// PRINT_INFO (option=0) no time info

#define PRINT_CONV (1)
// PRINT_INFO (option=1) time info

#define PRINT_PROC (1)

// #define PRINT_MED (1)
// print the field and interfaces passed between problems with MED interface


// ========= No user   ==================



// ============== no user ==================
#if PRINT_TIME==1
#include <ctime>
#endif



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


#endif // ------end file ---------------------