/*-------------------------------------------------------------------
 *
 ***  chemistryJointPdf    ***
 *
 * This file creates the joint PDF.
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/
#include "chemistryJointPdf.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
ChemistryJointPdf<V,M>::ChemistryJointPdf(const char* prefix, const VectorSet<V,M>& domainSet, int num_xi, int num_model_params, int num_catchalls)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"total").c_str(),domainSet), 
  n_xi(num_xi),
  n_k(num_model_params),
  n_atoms(num_catchalls)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54))
  {
    *m_env.subDisplayFile() << "Entering ChemistryJointPdf<V,M>::constructor()" << ": prefix = " << m_prefix << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54))
  {
    *m_env.subDisplayFile() << "Leaving ChemistryJointPdf<V,M>::constructor()" << ": prefix = " << m_prefix << std::endl;
  }
}

// Destructor --------------------------------------
template<class V,class M>
ChemistryJointPdf<V,M>::~ChemistryJointPdf()
{
}

// Math methods-------------------------------------
template<class V, class M>
double
ChemistryJointPdf<V,M>::actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  std::cout << "actual value is not implemented" << std::endl;
  if (domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal()) 
  {
    std::cout << "there is invalid input, this should be a queso error message" << std::endl;
  }

  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainDirection) {}; // just to remove compiler warning

  return 0.;
}
//--------------------------------------------------

template<class V, class M>
double
ChemistryJointPdf<V,M>::lnValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const
{
  if (gradVector   ) *gradVector     = m_domainSet.vectorSpace().zeroVector();
  if (hessianMatrix) *hessianMatrix *= 0.;
  if (hessianEffect) *hessianEffect  = m_domainSet.vectorSpace().zeroVector();

  if (domainVector[0]) {}; // just to remove compiler warning
  if (domainDirection) {}; // just to remove compiler warning

  double logpdf = 0.;
  double pdf = 0.;
  double pdfnorm = 0.;
  double pdflogn = 0.;

  // Mean of the prior mean for xi and the catchall reactions
  double meanMu = 3.0;  // mean
  double varMu = 0.01 * meanMu * meanMu;  // variance
  for (unsigned int i = 0; i < n_xi; i++)
  {
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0)) 
      { //  Set PDF to zero at \pm \infty
          pdf = 0.;
          logpdf = log(pdf);
      }
      else
      {
          pdfnorm = -log(sqrt(2 * M_PI * varMu)) - 1 / (2 * varMu) * pow((domainVector[i] - meanMu),2); 
          logpdf = logpdf + pdfnorm;
      }
  }

  // Variance of the prior variance for xi and the catchall reactions
  for (unsigned int i = n_xi; i < 2*n_xi; i++)
  {
      if ((domainVector[i] == -INFINITY ) || (domainVector[i] ==  INFINITY ) || (m_normalizationStyle != 0)) 
      {
         pdf = 0.;
         logpdf = log(pdf);
      }
      else
      {
          logpdf = logpdf - domainVector[i];
      }
  }

  // xi and the catchall reactions
  for (unsigned int i = 2*n_xi; i < 3*n_xi; i++)
  {
      double mu = domainVector[i - 2*n_xi];     // mean
      double eta = exp(domainVector[i - n_xi]); // variance
      double xi = exp(domainVector[i]);
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0))
      {
         pdf = 0.;
         logpdf = log(pdf);}
      else
      {
         pdflogn = -log(xi * sqrt(2 * M_PI * eta)) - 1 / (2 * eta) * pow((domainVector[i] - mu),2);
         logpdf = logpdf + pdflogn;
      }
  }
  // Catchall enthalpies:
  //   hk = \alpha_{0k} + \alpha_{1k}*T + \alpha_{2k} * T * T
  // Mean of the prior mean for \alpha_{0k}
  double meanMuEn = 10.0; // mean 
  double varMuEn = 0.01 * meanMuEn * meanMuEn;  // variance
  for (unsigned int i = 3*n_xi; i < 3*n_xi + n_atoms; i++)
  {
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0   ))
      {
         pdf = 0.;
         logpdf = log(pdf);
      }
      else
      {
         pdfnorm = -log(sqrt(2 * M_PI * varMuEn)) - 1 / (2 * varMuEn) *pow((domainVector[i] - meanMuEn),2); 
         logpdf = logpdf + pdfnorm;
      }
  }
  // Mean of the prior mean for \alpha_{ik}, i > 0
  // Note:  These are > 0 so we use a lognormal prior.
  for (unsigned int i = 3*n_xi + n_atoms; i < 3*n_xi + 3*n_atoms; i++)
  {
      double mu = -6.; //mean
      double eta = 0.01 * mu * mu;  //variance
      double beta = exp(domainVector[i]);
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0))
      {
         pdf = 0.;
         logpdf = log(pdf);}
      else
      {
         pdflogn = -log(beta * sqrt(2 * M_PI * eta)) - 1 / (2 * eta) * pow((domainVector[i] - mu),2); // with jacobian included
         logpdf = logpdf + pdflogn;
      }
  }
  // Variance of the prior variance for all catchall enthalpy coefficients
  for (unsigned int i = 3*n_xi + 3*n_atoms; i < 3*n_xi + 6*n_atoms; i++)
  {
      if ((domainVector[i] == -INFINITY ) || (domainVector[i] ==  INFINITY ) || (m_normalizationStyle != 0))
      {
         pdf = 0.;
         logpdf = log(pdf);
      }
      else
      {
         logpdf = logpdf - domainVector[i];
      }
  }
  // Catchall enthalpy coefficients for \alpha_{0k}
  for (unsigned int i = 3*n_xi + 6*n_atoms; i < 3*n_xi + 7*n_atoms; i++)
  {
      double mu = domainVector[i - 6*n_atoms]; //mean
      double eta = exp(domainVector[i - 3*n_atoms]); //variance
      double alpha = domainVector[i];
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0))
      {
         pdf = 0.;
         logpdf = log(pdf);}
      else
      {
         pdflogn = -log(sqrt(2 * M_PI * eta)) - 1 / (2 * eta) * pow((alpha - mu),2); 
         logpdf = logpdf + pdflogn;
      }
  }
  // Catchall enthalpy coefficients for \alpha_{ik}, i > 0
  for (unsigned int i = 3*n_xi + 7*n_atoms; i < 3*n_xi + 9*n_atoms; i++)
  {
      double mu = domainVector[i - 6*n_atoms]; // mean
      double eta = exp(domainVector[i - 3*n_atoms]); // variance
      double beta = exp(domainVector[i]);
      if ((domainVector[i] == -INFINITY) || (domainVector[i] ==  INFINITY) || (m_normalizationStyle != 0))
      {
         pdf = 0.;
         logpdf = log(pdf);}
      else 
      {
         pdflogn = -log(beta * sqrt(2 * M_PI * eta)) - 1 / (2 * eta) * pow((domainVector[i] - mu),2); // with jacobian included
         logpdf = logpdf + pdflogn;
      }
  }
/*
  // Global activation energy
  double mu_Tg = 1000.0; // Mean
  double var_Tg = 0.01 * mu_Tg * mu_Tg; // Variance
  pdfnorm = -log(sqrt(2.0 * M_PI * var_Tg)) - 1.0 / (2.0 * var_Tg) * pow((domainVector[3*n_xi + n_k + 9*n_atoms] - mu_Tg), 2.0);
  logpdf  = logpdf + pdfnorm;
*/

  double mu_Tigg = 919.0; // Mean
  double var_Tigg = 0.01 * mu_Tigg * mu_Tigg; // Variance
  pdfnorm = -log(sqrt(2.0 * M_PI * var_Tigg)) - 1.0 / (2.0 * var_Tigg) * pow((domainVector[3*n_xi + 9*n_atoms] - mu_Tigg), 2.0);
  logpdf  = logpdf + pdfnorm;

  double mu_Tadg = 2680.0; // Mean
  double var_Tadg = 0.01 * mu_Tadg * mu_Tadg; // Variance
  pdfnorm = -log(sqrt(2.0 * M_PI * var_Tadg)) - 1.0 / (2.0 * var_Tadg) * pow((domainVector[3*n_xi + 9*n_atoms + 1] - mu_Tadg), 2.0);
  logpdf  = logpdf + pdfnorm;

  // Now return the joing pdf
  return logpdf;
}
//--------------------------------------------------
template<class V, class M>
double
ChemistryJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) 
  {
    *m_env.subDisplayFile() << "Entering ChemistryJointPdf<V,M>::computeLogOfNormalizationFactor()" << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2))
  {
    *m_env.subDisplayFile() << "Leaving ChemistryJointPdf<V,M>::computeLogOfNormalizationFactor()" << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54))
  {
    *m_env.subDisplayFile() << " return normalization factor " << std::endl;
  }
  return value;
}

}  // End namespace QUESO

template class QUESO::ChemistryJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
