// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XFIT_H_
#define _AUROSTD_XFIT_H_

//********************************************************************************
//              Functions to work with polynomials
namespace aurostd{
  double evalPolynomial(double x, const xvector<double> &p);
  void evalPolynomialDeriv(double x, const xvector<double> &p, xvector<double> &dp);
  xvector<double> evalPolynomialDeriv(double x, const xvector<double> &p, uint n);
  xmatrix<double> Vandermonde_matrix(const xvector<double> &x, int n);
  double polynomialFindExtremum(const xvector<double> &p, double xmin, double xmax,
      double tol=_mm_epsilon);
}

//********************************************************************************
//              Definitions of the NonlinearFit class members
namespace aurostd{
  /// Fit to a nonlinear model using the Levenberg-Marquardt algorithm.
  ///
  /// The implementation here is based on the ideas from Numerical Recipes and
  /// K. Madsen et al. Methods For Non-linear Least Squares Problems
  /// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  ///
  /// Caution: the default value for the parameter tau (a scaling factor for the initial step size)
  /// was picked to yield a correct fit to the Murnaghan equation of state. 
  /// If you observe that this is not a good choice for your function,
  /// try 1e-6 if the initial guess is believed to be a good approximation to the true
  /// parameters. Otherwise 1e-3 or even 1 might be a better choice.
  class NonlinearFit{
    public:
      NonlinearFit();
      NonlinearFit(const NonlinearFit &nlf);
      NonlinearFit(xvector<double> &x, xvector<double> &y, xvector<double> &guess,
          double foo(const double x, const xvector<double> &p, xvector<double> &dydp),
          double tol=1e-6, double tau=1e-12, int max_iter=1000);
      ~NonlinearFit();
      const NonlinearFit& operator=(const NonlinearFit &qha);
      int Npoints, Nparams;
      double tol; /// convergence tolerance criterion
      double tau; /// scaling parameter for initial step size
      int max_iter; /// maximum number of allowed iterations
      xvector<double> x,y;   // data points
      xvector<double> residuals; // residuals of a given model function
      xvector<double> guess; // initial guess for fit parameters
      xvector<double> p;     // parameters obtained by fit
      xvector<double> dydp;  // derivative of a given function w.r.t parameters
      xmatrix<double> A;     // J^T.J matrix
      xmatrix<double> J;     // Jacobian of a model function w.r.t parameters
      double (*f)(const double x, const xvector<double> &p, xvector<double> &dydp);
      bool fitLevenbergMarquardt();
      void Jacobian(const xvector<double> &guess);
      void calculateResiduals(const xvector<double> &params);
      double calculateResidualSquareSum(const xvector<double> &params);
      void clear();

    private:
      void free();
      void copy(const NonlinearFit &nlf);
  };
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

