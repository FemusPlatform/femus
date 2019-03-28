The turbulence library is made of three different blocks:

- Utility 
- Turbulence model solvers
- Systems for \f$\nu_t\f$ and \f$\alpha_t\f$

The interaction between the different blocks is represented below

\dot

digraph D {
    compound=true;
    subgraph cluster_21 {
      label="TurbUtils";
      
      
      subgraph cluster_dyn {
        label="_DynModel"; 
        URL="\ref #TurbUtils::_DynModel";
        
    
        subgraph cluster_4pdin {
          label="4 par. model"; 
          URL="\ref #Dyn4P";
          style=filled;
          color=lightseagreen;
    
           a11 [ label=<k &#949;>, URL="\ref #Dyn4P_KE"];
           a12 [ label=<k &#969;>, URL="\ref #Dyn4P_KW"];
           a13 [ label=<k &#937;>, URL="\ref #Dyn4P_logKW"];
        }
        
        subgraph cluster_wilcox {
          label = "Wilcox model"; 
          URL="\ref #DynWilcox";
          style=filled;
          color=lightseagreen;
    
           b11 [ label=<k &#969;>, URL="\ref #DynWilcox_KW"];
           b12 [ label=<k &#937;>, URL="\ref #DynWilcox_logKW"];    
        }
      }
      
      subgraph cluster_therm {
        label="_ThermModel"; 
        URL="\ref #TurbUtils::_ThermModel";
        
    
        subgraph cluster_4ptherm {
          label="4 par. model"; 
          URL="\ref #Therm4P";
          style=filled;
          color=red;
    
           c11 [ label=<k<sub>&#952;</sub> &#949;<sub>&#952;</sub>>, URL="\ref #Therm4P_KHEH"];
           c12 [ label=<k<sub>&#952;</sub> &#969;<sub>&#952;</sub>>, URL="\ref #Therm4P_KHWH"];
           c13 [ label=<k<sub>&#952;</sub> &#937;<sub>&#952;</sub>>, URL="\ref #Therm4P_logKHWH"];
        }
      }
    }
    
    subgraph cluster_2 {
        label=<Solvers for &#957;<sub>t</sub> and &#945;<sub>t</sub>>;
        e [ label=<&#957;<sub>t</sub>> URL="\ref #MGSolTURB"];
        f [ label=<&#945;<sub>t</sub>> URL="\ref #MGSolTURB"];
    }
    subgraph cluster_1 {
        label="Turbulence model solvers";
       c [ label="Dynamical" URL="\ref #MGSolRANS"];
       d [ label="Thermal" URL="\ref #MGSolRANS_thermal"];
    }
    
    a12 -> e [ltail=cluster_21, lhead=cluster_1, label="CalcMuT", URL="\ref #DYNturModels::CalcMuT" ]; 
    c12 -> f [ltail=cluster_therm, lhead=f, label="CalcAlphaT", URL="\ref #Therm4P::CalcAlphaTurb" ];
    a12 -> c [ltail=cluster_21, lhead=cluster_1, label="Source terms and near wall values"]; 
}

\enddot

#### Utility blocks

The utility block is used to compute the eddy diffusion coefficients \nu_t and \alpha_t,
the turbulence models source and dissipative terms and other variables of interest, as 
friction velocity u_\f$\tau\f$.
The idea is that the turbulence model solvers interact with the utility block for the calculation
of the terms needed to solve the model equations. This helps reducing redundancy of implemented
functions and a more easy control of modeled terms.
Turbulence model solvers interact with the utility block through a #TurbUtils class object.
All the solvers interact with the *same* #TurbUtils class object by using #MGUtils::_TurbParameters,
which is a pointer to a #TurbUtils class object that is built at main level.
For each implemented turbulence model a structure for the calculation of turbulence variables
near wall behaviors, turbulence variable initial values, eddy diffusivity coefficients, source and
dissipative terms is provided.
They are organized as follows

- #DYNturModels: Top level structure for dynamical turbulence models
  - #Dyn4P: Specialization for the four parameter turbulence model
     - #Dyn4P_KE: Specialization for the \f$k ~\varepsilon\f$ version of the four parameter model
     - #Dyn4P_KW: Specialization for the \f$k ~\omega\f$ version of the four parameter model
     - #Dyn4P_logKW: Specialization for the logarithmic \f$K ~\Omega\f$ version of the four parameter model
  - #DynWilcox: Specialization for the  \f$k ~\omega\f$ Wilcox turbulence model
     - #DynWilcox_KW: Original \f$k ~\omega\f$ Wilcox turbulence model
     - #DynWilcox_logKW: Logarithmic formulation of the \f$k ~\omega\f$ Wilcox turbulence model
- #Therm4P: Top level structure for thermal turbulence models, based on four parameter model
  - #Therm4P_KHEH: Specialization for the \f$k_\theta ~\varepsilon_\theta\f$ version of the four parameter model
  - #Therm4P_KHWH: Specialization for the \f$k_\theta ~\omega_\theta\f$ version of the four parameter model
  - #Therm4P_logKHWH: Specialization for the logarithmic \f$K_\theta ~\Omega_\theta\f$ version of the four parameter model

The TurbUtils class is a wrapped around #DYNturModels, #Therm4P and their derived structure 
functions. A pointer to both #DYNturModels and #Therm4P classes is contained inside TurbUtils
class, as #TurbUtils::_DynModel and #TurbUtils::_ThermModel respectively. The proper structures
are built during code execution, depending on the selected turbulence models.

#### Systems for \f$\nu_t\f$ and \f$\alpha_t\f$



