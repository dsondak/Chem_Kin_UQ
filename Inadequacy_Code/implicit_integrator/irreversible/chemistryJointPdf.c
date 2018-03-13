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
#include <queso/VectorSpace.h>
#include <queso/BasicPdfsBase.h>

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
  queso_require_equal_to_msg(domainVector.sizeLocal(), this->m_domainSet.vectorSpace().dimLocal(), "invalid input");

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
  std::vector<double> 
    hypermean_mean = {log(2.2716495127060917), 5.2756035602390422, 1.1021099751123025, 
                      log(3.7714830544483780), 1.0e-16, 8.9433281629646342, 
                      log(1.7721512823487061), 3.8246742586538279, 2.0745500053927372, 
                      log(5.5706916507985680), 1.0e-16, 1.3135622105611983, 
                      log(8.2830919705299377), -1.3585766873971914, 1.3005166005172624, 
                      log(3.5305061634938347), 1.0e-16, 1.7168401051086030,
                      log(6.6857188699084237), 6.6634385904742177, 7.9705662739958876, 
                      log(7.8430783348457687), 1.0e-16, 1.0228335042104413,
                      log(4.3298112708110565), -6.9922795750273214, 9.0747132030759138, 
                      log(3.7484557650026878), 1.0e-16, 4.3359623996979244};

  double hypermean_variance; // Variance of the hypermean for inadequacy kinetics
  double var = 0.01; // Variance factor (so std = 10% of mean)
  double var1 = 0.5; // Variance factor for modified Arrhenius parameters
  int p_hypermean_kinetics = 3 * n_reactions_inad; // Pointer to index
  for (unsigned int i = 0; i < 3 * n_reactions_inad; i++)
  {
      if (hypermean_mean[i] != 1.0e-16)
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

  // Variance of the prior variance for inadequacy parameters
  const unsigned int n_params = 3 * n_reactions_inad;
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

  printf("\n  DomainVector \n");
  for (unsigned int i = 0; i < 9 * n_reactions_inad; i++) { 
      printf("%25.16e\n", domainVector[i]);
  }
  printf("\n");
  std::cout << "\n joint prior = " << logpdf << "\n" << std::endl;
  //if (logpdf <= -1.0e+03) {
  //   exit(0);
  //}

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
