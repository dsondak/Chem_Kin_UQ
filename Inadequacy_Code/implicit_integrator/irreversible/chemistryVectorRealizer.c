#include "chemistryVectorRealizer.h"
#include <limits>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

//only returning vector of one's because it wouldn't let me return NULL
//shouldn't be actually used in chemistry code
namespace QUESO {

// Constructor -------------------------------------
template<class V, class M>
ChemistryVectorRealizer<V,M>::ChemistryVectorRealizer(
  const char*                  prefix,
  const VectorSet<V,M>& unifiedImageSet)
  :
  BaseVectorRealizer<V,M>(((std::string)(prefix)+"gen").c_str(),unifiedImageSet,std::numeric_limits<unsigned int>::max())
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering ChemistryVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving ChemistryVectorRealizer<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V, class M>
ChemistryVectorRealizer<V,M>::~ChemistryVectorRealizer()
{
}
// Realization-related methods----------------------
template<class V, class M>
void
ChemistryVectorRealizer<V,M>::realization(V& nextValues) const
{
 /* 
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);

  if (imageBox == NULL) {
    queso_error_msg(" For ChemistryVectorRealizer<V,M>::realization(), only box images are supported right now ");}
  //take log of Chemistry bounds to set uniform bounds
  GslVector logMinValues(imageBox->minValues());
  for (unsigned int i = 0; i < logMinValues.sizeLocal(); ++i) {
    if (logMinValues[i] < 0.) {
      queso_error_msg("The minimum value for a Chemistry distribution should be greater than or equal to zero."); }
    else logMinValues[i] = std::log(logMinValues[i]);
  }

  GslVector logMaxValues(imageBox->maxValues());
  for (unsigned int i = 0; i < logMaxValues.sizeLocal(); ++i) {
    if (logMaxValues[i] <= 0.) {
      queso_error_msg("The maximum value for a Chemistry distribution should be greater than zero."); }
    else logMaxValues[i] = std::log(logMaxValues[i]);
  }


  nextValues.cwSetUniform(logMinValues,logMaxValues);
  for (unsigned int i = 0; i < nextValues.sizeLocal(); ++i){
    nextValues[i] = exp(nextValues[i]);
  }
 */
  const BoxSubset<V,M>* imageBox = dynamic_cast<const BoxSubset<V,M>* >(&m_unifiedImageSet);
  GslVector minValues(imageBox->minValues());
  minValues.cwSet(-5.);
  GslVector maxValues(imageBox->maxValues());
  maxValues.cwSet(5.);
  nextValues.cwSetUniform(minValues,maxValues);
  //for (unsigned int i = 0; i < nextValues.sizeLocal(); i++){
    //nextValues[i].setUniform(-5.,5.);
  //}
  //nextValues.cwSet(1.0);
  return;
}

}  // End namespace QUESO

template class QUESO::ChemistryVectorRealizer<QUESO::GslVector, QUESO::GslMatrix>;
