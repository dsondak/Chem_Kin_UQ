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
ChemistryJointPdf<V,M>::ChemistryJointPdf(const char* prefix, const VectorSet<V,M>& domainSet, 
                                          int num_reactions_inad, int num_atoms)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"total").c_str(),domainSet), 
  n_reactions_inad(num_reactions_inad),
  n_atoms(num_atoms)
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

  double pdf    = 0.0;
  double logpdf = 0.0;

  // Mean of the prior mean for inadequacy kinetics
  std::vector<double> hypermean_mean = {1.0, 0.0, 1.65,
                                        2.0, 0.0, 1.65, 
                                        1.0, 0.0, 1.65,
                                        1.5, 0.0, 1.65, 
                                        1.17, 0.0, 1.65};
  double hypermean_variance; // Variance of the hypermean for inadequacy kinetics
  double var = 0.01; // Variance factor (so std = 10% of mean)
  double var1 = 0.5; // Variance factor for modified Arrhenius parameters
  int p_hypermean_kinetics = 3 * n_reactions_inad + 4 * n_atoms; // Pointer to index
  for (unsigned int i = 0; i < 3 * n_reactions_inad; i++)
  {
      if (hypermean_mean[i] != 0)
      {
         hypermean_variance = var * hypermean_mean[i] * hypermean_mean[i];
      }
      else
      {
         hypermean_variance = var1; // Only for modified Arrhenius parameters
      }
      pdf = -log(sqrt(2.0 * M_PI * hypermean_variance)) - 
                 1.0 / 2.0 / hypermean_variance * pow((domainVector[p_hypermean_kinetics + i] - hypermean_mean[i]), 2); 
      logpdf += pdf;
  }

  // Mean of the prior mean for inadequacy thermochemistry
  hypermean_mean = {3.90372558e+04, 13.6559654, 1.20459536e-03, 5.0, 
                    3.90372558e+04, 13.6559654, 1.20459536e-03, 10.0};
  p_hypermean_thermo = 6 * n_reactions_inad + 4 * n_atoms; // Pointer to index
  for (unsigned int i = 0; i < 4 * n_atoms; i++)
  {
      if (hypermean_mean[i] != 0)
      {
         hypermean_variance = var * hypermean_mean[i] * hypermean_mean[i];
      }
      else
      {
         hypermean_variance = var;
      }
      pdf = -log(sqrt(2.0 * M_PI * hypermean_variance)) - 
                 1.0 / 2.0 / hypermean_variance * pow((domainVector[p_hypermean_thermo + i] - hypermean_mean[i]), 2); 
      logpdf += pdf;
  }

  // Variance of the prior variance for inadequacy parameters
  const unsigned int n_params = 3 * n_reactions_inad + 4 * n_atoms;
  int p_hypervar = 2 * n_params;
  for (unsigned int i = 0; i < n_params; i++)
  {
      logpdf -= domainVector[p_hypervar + i]; // Jeffries prior
  }

  // Inadequacy kinetics parameters
  double mean;     // mean
  double variance; // variance
  double theta;    // parameter value
  for (unsigned int i = 0; i < 3 * n_reactions_inad; i++)
  {
      mean     = domainVector[p_hypermean_kinetics + i];
      variance = exp(domainVector[p_hypervar + i]);
      if (i % 3 == 0)
      {
         theta = exp(domainVector[i]);
         pdf   = -log(theta * sqrt(2.0 * M_PI * variance)) - 1.0 / 2.0 / variance * pow((domainVector[i] - mean), 2);
      }
      else
      {
         pdf = -log(sqrt(2.0 * M_PI * variance)) - 1.0 / 2.0 / variance * pow((domainVector[i] - mean), 2);
      }

      logpdf += pdf;
  }

  // Inadequacy thermo chemistry parameters
  double mean;     // mean
  double variance; // variance
  double theta;    // parameter value
  for (unsigned int i = 0; i < 4 * n_atoms; i++)
  {
      mean     = domainVector[p_hypermean_thermo + i];
      variance = exp(domainVector[p_hypervar + 3 * n_reactions_inad + i]);
      if (i % 4 == 0)
      {
         pdf  = -log(sqrt(2.0 * M_PI * variance)) - 1.0 / 2.0 / variance * pow((domainVector[3 * n_reactions_inad + i] - mean), 2);
      }
      else
      {
         theta = exp(domainVector[3 * n_reactions_inad + i]);
         pdf   = -log(theta * sqrt(2.0 * M_PI * variance)) - 1.0 / 2.0 / variance * pow((domainVector[3 * n_reactions_inad + i] - mean), 2);
      }

      logpdf += pdf;
  }

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
