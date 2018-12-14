#ifndef __EquationMapFiller_h__
#define __EquationMapFiller_h__

#include <vector>

class FEMUS;
class EquationSystemsExtendedM;

enum FIELDS : int;

class EquationMapFiller{  
    
public:    
    int _nvars[3];
    
    EquationMapFiller(){}
    ~EquationMapFiller(){}
    void FillEquationMap(FEMUS & FemusProblem, const std::vector<FIELDS> & pbName);
    
    void setNavierStokes(EquationSystemsExtendedM & EqMap);
    void setTemperature(EquationSystemsExtendedM & EqMap);
    void setDynamicTurbulence(EquationSystemsExtendedM & EqMap);
    void setThermalTurbulence(EquationSystemsExtendedM & EqMap);
    void setFluidStructure(EquationSystemsExtendedM & EqMap);
    void setColor(EquationSystemsExtendedM & EqMap);
    void setDisplacements(EquationSystemsExtendedM & EqMap);
    void setImmersedBoundary(EquationSystemsExtendedM & EqMap);
};

#endif
