#ifndef __Utils_DynTurModels_h__
#define __Utils_DynTurModels_h__

#include <math.h>
#include <algorithm>


using namespace std;

/**
 * Father structure for eddy diffusivity (\f$\nu_t\f$ and \f$\alpha_t \f$) coefficient calculation.
 */
struct DYNturModels {
    /// Molecular kinematic viscosity
    double _nu;

    /// Flag for Park limit
    int _Park;
    /// Flag for Yap source term
    int _YapCorr;

    /// Vector used to store real \f$k\f$ and \f$\omega \f$ variables
    double _KappaAndOmega[2];
    
    double _LowerKappa, _LowerOmega;
    
    DYNturModels () :_Park ( 0 ), _YapCorr ( 0 ), _LowerKappa (1.e-7), _LowerOmega (1.e-7)
    {
    };

    /// Function used to calculate eddy kinematic viscosity \f$ \nu_t \f$.
    /**
     *  The function is specialized by #Dyn4P and #DynWilcox derived structures
     */
    virtual double CalcMuT (
        double *TurbVariables, /**< Model dynamic turbulence variables */
        double ydist = 1.e-5   /**< Wall distance */
    )
    {
    };

    /// Function used to calculate source and siddipation terms for turbulence model equations
    /**
     *  The function is specialized by #Dyn4P and #DynWilcox derived structures, and their
     *  further derived structures for real or logarithmic turbulence variables
     */
    virtual void CalcDynSourceTerms (
        double *TurbVariables, /**< Model dynamic turbulence variables*/
        double *SourceTerms,   /**< Array used to store source terms  */
        double *DissTerms,     /**< Array used to store dissipative terms*/
        double MuT,            /**< Eddy kinematic viscosity, input value from #MGSolTURB solver */
        double VelGradMod,     /**< Modulus of velocity gradient */
        double ydist = 1.e-5   /**< Wall distance */
    )
    {
    };

    /// Function used to set molecular kinematic viscosity \f$\nu\f$ value
    void inline SetNu ( double nu )
    {
        _nu = nu;
    };

    /// Function used to set the value of Park limit flag
    void inline SetPark ( int Park )
    {
        _Park = Park;
    }

    /// Function used to set the flag value for Yap source term calculation
    void inline SetYap ( int Yap )
    {
        _YapCorr = Yap;
    }

    /// Function used to calculate Yap source term
    /**
     *  The function returns the Yap source term for \f$\varepsilon\f$ equation. Derived
     *  structures then adapt it to the specific turbulence model
     */
    double YapTerm (
        double *KappaAndEpsilon, /**< Real \f$k\f$ and \f$\varepsilon\f$ turbulence variables */
        double ydist = 1.e-5     /**< Wall distance */
    );

    /// Function used to calculate limited production terms for real \f$k\f$ equation
    double inline Park_KSource (
        double prod_k,           /**< \f$k\f$ original source term */
        double kappa,            /**< \f$k\f$ value */
        double VelGradMod        /**< Velocity gradient tensor modulus */
    )
    {
        return min ( prod_k, sqrt ( 8. / 3. ) * kappa * sqrt ( VelGradMod ) );
    };

    /// Function used to calculate near wall behaviors of real \f$k\f$ and \f$\omega\f$ variables
    void DynTurNearKappaAndOmegaValues (
        double *TurbValues,     /**< Input array used to store real \f$k\f$ and \f$\omega\f$ near wall values */
        double WallDist,        /**< Wall distance */
        double Utau             /**< Friction velocity \f$u_\tau\f$ */
    );

    /// Function used to convert real \f$k\f$ and \f$\omega\f$ values to specific turbulence model variable values
    /**
     *  Real \f$k\f$ and \f$\omega\f$ are read from #DYNturModels::_KappaAndOmega array
     */
    virtual void ConvertKandWtoLocal (
        double *TurbValues      /**< Input array used to store converted variable values */
    ) {};

    /// Function used to compute real \f$k\f$ and \f$\omega\f$ values from specific turbulence model variable values
    /**
     *  Computed values are stored inside #DYNturModels::_KappaAndOmega array
     */
    virtual void CalcKappaAndOmega (
        double *TurbVariables   /**< Model dynamical turbulence variables */
    ) {};

    /// Function used to calculate initial values of real \f$k\f$ and \f$\omega\f$ variables
    void DynTurInitKappaAndOmegaValues (
        double *TurbValues,     /**< Input array used to store converted variable values */
        double vel,             /**< Reference velocity modulus */
        double diameter         /**< Reference diameter */
    );

};

#endif

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
