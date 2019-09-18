#ifndef __mgtimeloop_h__
#define __mgtimeloop_h__

// Forward class
class MGUtils;
class EquationSystemsExtendedM;

// ===============================================
//                  MGTimeLoop class
// ===============================================
/// Class doing generic time loop for the solution of the system
class MGTimeLoop {
 protected:
  // Data ---------------------------
  MGUtils& _mgutils;                   ///< MGUtils vector of pointers
  EquationSystemsExtendedM& _mgeqmap;  ///< MGEquationsMap vector of pointers

 public:
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  /// This function is the constructor
  MGTimeLoop(
      MGUtils& _mgutils,                  ///< MGUtils class (data)
      EquationSystemsExtendedM& _mgeqmap  ///< EquationSystemsExtendedM
  );
  /// This function is the class destructor
  ~MGTimeLoop(){};
  ///@}

  //-----------------------------------------------------------------------------
  ///@{ \name FUNCTIONS FOR TIME LOOP
  /// This function setup a time loop (t=0)
  void transient_setup(
      const int& restart,  ///< restart iteration (0=no restart)  (in)
      int& t_in,           ///< initial time iteration            (in)
      double& time         ///< running time                      (in)
  );
  // ========================================================================================
  /// This function runs a time loop
  void transient_loop(
      const int& t_in,    ///< initial time iteration          (in)
      double& time,       ///< running time                    (in)
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );

  // ========================================================================================
  /// This function runs a time loop
  void transient_loop_no_up(
      const int& t_in,    ///< initial time iteration          (in)
      double& time,       ///< running time                    (in)
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );

  void transient_update(
      const int& t_in,    ///< initial time iteration          (in)
      double& time,       ///< running time                    (in)
      const int& eq_min,  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max   ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );

  // ========================================================================================
  /// This function runs  a time step
  void transient_onestep(
      const int& t_in,        ///< initial time iteration      (in)
      const int& t_step,      ///< running time iteration      (in)
      const int& print_step,  ///< print every                 (in)
      double& time,           ///< running time                (in)
      double& dt,             ///< step time                   (in)
      const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max       ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );
  // ========================================================================================
  /// This function runs  a time step
  void transient_solve_and_update(
      const int& t_in,        ///< initial time iteration      (in)
      const int& t_step,      ///< running time iteration      (in)
      const int& print_step,  ///< print every                 (in)
      double& time,           ///< running time                (in)
      double& dt,             ///< step time                   (in)
      const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max       ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );
  // ========================================================================================
  /// This function runs  a time step
  void steady(
      const int& nmax_step,   ///< number max of steps
      const double& toll,     ///< tolerance
      const int& t_step,      ///< running time iteration      (in)
      const int& print_step,  ///< print every                 (in)
      const double& dt,       ///< step time                   (in)
      const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h) (in)
      const int& eq_max       ///< eq max to solve -> enum  FIELDS (equations_conf.h) (in)
  );

  void set_uooold(
      const int& vec_from,  ///< source vector to be copied      (in)
      const int& vec_to,    ///< target vector                   (in)
      const double& toll,   ///< running time                    (in)
      const double& time,   ///< running time                    (in)
      const int& eq_min,    ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max     ///< eq max to solve -> enum  FIELDS (equations_conf.h)
  );

  // ========================================================================================
  /// This function does one step of control solution
  void transient_control_onestep(
      const int& nmax_step,   ///< number max of steps         (in)
      const int& it,          ///< iteration                   (in)
      const int& t_step,      ///< running time iteration      (in)
      const int& print_step,  ///< print every                 (in)
      double& time,           ///< running time                (in)
      double& dt,             ///< step time                   (in)
      const int& eq_min,      ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max,      ///< eq max to solve -> enum  FIELDS (equations_conf.h)
      bool& converged         ///< check if the solution converged (1->converged)     (out)
  );
  // ========================================================================================
  // ========================================================================================
  /// This function does one step of control solution
  void transient_control_onestep(
      const int& nmax_step,               ///< number max of steps         (in)
      const int& it,                      ///< iteration                   (in)
      const int& t_step,                  ///< running time iteration      (in)
      const int& print_step,              ///< print every                 (in)
      double& time,                       ///< running time                (in)
      double& dt,                         ///< step time                   (in)
      const int& eq_min,                  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max,                  ///< eq max to solve -> enum  FIELDS (equations_conf.h)
      std::vector<double> controlled_eq,  /// vector with the number of convergence-controlled equation
      bool& converged,                    ///< check if the solution converged (1->converged)     (out)
      const double& toll);
  // ========================================================================================
  /// This function runs a time step with a feedback on the norm convergence
  void transient_underrelaxed_onestep(
      const int& it,                      ///< initial time iteration      (in)
      const int& t_step,                  ///< running time iteration      (in)
      const int& print_step,              ///< print every                 (in)
      double& time,                       ///< running time                (in)
      double& dt,                         ///< step time                   (in)
      const int& eq_min,                  ///< eq min to solve -> enum  FIELDS (equations_conf.h)
      const int& eq_max,                  ///< eq max to solve -> enum  FIELDS (equations_conf.h)
      std::vector<double> controlled_eq,  /// vector with the number of convergence-controlled equation
      bool& converged,                    ///< check if the solution converged (1->converged) (out)
      const double& toll);
  // ========================================================================================
  /// This function  prints a file time.****.xmf in RESU
  void transient_print_xmf(int t_init  ///< initial time iteration     (in)
  );
  ///@}
  // ===========================================================================
  /// This function controls the transient loop
  void dummy_step(
      const int& t_in,        ///< initial time iteration      (in)
      const int& t_step,      ///< running time iteration      (in)
      const int& print_step,  ///< print every                 (in)
      double& time,           ///< running time                (in)
      double& dt              ///< step time                   (in)
  );
};

#endif  //  ---------  end header -------------------------
