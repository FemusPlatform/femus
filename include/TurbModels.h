#ifndef __TURB_MODELS__
#define __TURB_MODELS__

#include <string>
#include <map>
  enum DynTurbModel{
    // NAGANO BASED MODELS
    nagano_ke   =1,
    nagano_kw   =2,
    nagano_log  =3,
    wilcox      =4,
    wilcox_log  =5,
    default_dyn =0,
    nagano_k    =10,
    nagano_logk =11,
    nagano_w    =10,
    nagano_logw =11,
    nagano_e    =12,
    // WILCOX BASED MODELS
    wilcox_k    =20,
    wilcox_logk =21,
    wilcox_nut  =22,
    wilcox_w    =20,
    wilcox_logw =21
  };
  
  
  enum ThermTurbModel{
    default_therm = 0,
    nagano_keT = 1,
    nagano_kwT = 2,
    nagano_logT = 3,
    natural_kh = 10,
    logarithmic_kh = 11,
    natural_omegah = 20,
    logarithmic_omegah = 21,
    kays = 30,
  };


  const std::map<std::string, DynTurbModel> DynTurbModelMap={
    {"nagano_ke"  ,     nagano_ke  },
    {"nagano_kw"  ,     nagano_kw  },
    {"nagano_log" ,     nagano_log },
    {"wilcox"     ,     wilcox     },
    {"wilcox_log" ,     wilcox_log },
    {"default"    ,     default_dyn}
  };
  const std::map<std::string, ThermTurbModel> ThermTurbModelMap={
    {"nagano_keT"  ,     nagano_keT  },
    {"nagano_kwT"  ,     nagano_kwT  },
    {"nagano_logT" ,     nagano_logT },
    {"kays"        ,     kays        },
    {"default"     ,     default_therm}
  };

  
#endif
