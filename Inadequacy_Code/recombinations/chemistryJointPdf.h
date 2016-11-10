#ifndef UQ_CHEMISTRY_JOINT_PROB_DENSITY_H
#define UQ_CHEMISTRY_JOINT_PROB_DENSITY_H

#include <queso/JointPdf.h>
#include <queso/UniformJointPdf.h> //actually not using this or following line yet, but should
#include <queso/JeffreysJointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

/*!
 * \class ChemistryJointPdf
 * \brief A class for handling chemistry joint PDFs.
 *
 * This class allows the mathematical definition of a Chemistry Joint PDF.
*/

template<class V, class M>
class ChemistryJointPdf : public BaseJointPdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object of the class, given a prefix and the domain set of the jeffreys PDF.  */
  ChemistryJointPdf(const char* prefix,
                    const VectorSet<V,M>& totalDomain,
                    int num_reactions_inad, int num_atoms);
  //! Destructor
 ~ChemistryJointPdf();
 //@}

   //! @name Math methods
  //@{
  //! Actual value of the jeffreys PDF.
  /*! If the domain of the PDF is well defined (neither negative nor infinite), then the actual
   * value is given by 1/x, otherwise the actual
   * value is 0.*/
  
  double actualValue(const V& domainVector, const V* domainDirection, V*
      gradVector, M* hessianMatrix, V* hessianEffect) const;
  
  //! Logarithm of the value of the jeffreys PDF.
  /*! Analogous to the actualValue routine, except that the logarithm of the calculated value is
   * returned. */
  double lnValue (const V& domainVector, const V* domainDirection, 
                  V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  int n_reactions_inad;
  int n_atoms;
  
  //TODO: do we want this part?
  //! Computes the logarithm of the normalization factor.
  /*! This routine calls BaseJointPdf::commonComputeLogOfNormalizationFactor().*/
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const;
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_normalizationStyle;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;
};

}  // End namespace QUESO

#endif // UQ_CHEMISTRY_JOINT_PROB_DENSITY_H
