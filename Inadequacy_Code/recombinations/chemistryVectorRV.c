#include "chemistryVectorRV.h"
#include "chemistryJointPdf.h"
#include "chemistryVectorRealizer.h"
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor-------------------------------
template<class V, class M>
ChemistryVectorRV<V,M>::ChemistryVectorRV(const char* prefix, const VectorSet<V,M>& imageSet, int num_xi, int num_model_params, int num_catchalls)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"total").c_str(),imageSet),
  n_xi(num_xi),
  n_k(num_model_params),
  n_atoms(num_catchalls)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54))
  {
    *m_env.subDisplayFile() << "Entering ChemistryVectorRV<V,M>::constructor()" << ": prefix = " << m_prefix << std::endl;
  }
  
  // the pdf is the only thing we need
  m_pdf        = new ChemistryJointPdf<V,M>(m_prefix.c_str(), m_imageSet, n_xi, n_k, n_atoms);
  m_realizer   = NULL; // new ChemistryVectorRealizer<V,M>(m_prefix.c_str(),m_imageSet);
  m_subCdf     = NULL; // FIX ME: complete code
  m_unifiedCdf = NULL; // FIX ME: complete code
  m_mdf        = NULL; // FIX ME: complete code

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54))
  {
    *m_env.subDisplayFile() << "Leaving ChemistryVectorRV<V,M>::constructor()" << ": prefix = " << m_prefix << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
ChemistryVectorRV<V,M>::~ChemistryVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods --------------------------------------
template <class V, class M>
void
ChemistryVectorRV<V,M>::print(std::ostream& os) const
{
  os << "ChemistryVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::ChemistryVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
