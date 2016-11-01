#ifndef UQ_CHEMISTRY_VECTOR_RV_H
#define UQ_CHEMISTRY_VECTOR_RV_H

#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/JointPdf.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorCdf.h>
#include <queso/VectorMdf.h>
#include <queso/SequenceOfVectors.h>
#include <gsl/gsl_sf_psi.h> // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

namespace QUESO {

/*!
 * \class ChemistryVectorRV
 * \brief A class representing a chemistry vector RV.
 * 
 * This class allows the user to compute the value of my custom chemistry PDF. 
 * Right now, won't worry about generating realizations
 * (samples) from it. */

template<class V, class M>
class ChemistryVectorRV : public BaseVectorRV<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor
  /*! Constructs a uniform vector RV, given a prefix and the image set of the vector RV.*/
  ChemistryVectorRV(const char*                  prefix,
//                         const VectorSet<V,M>& imageSetMu,
//                         const VectorSet<V,M>& imageSetSigma,
//                         const VectorSet<V,M>& imageSetC,
//                         const VectorSet<V,M>& imageSetK,
                         const VectorSet<V,M>& imageSet,
                         int num_xi, int num_model_params, int num_catchalls);
  //! Virtual destructor
  virtual ~ChemistryVectorRV();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

  int n_xi;
  int n_k;
  int n_atoms;

private:
  using BaseVectorRV<V,M>::m_env;
  using BaseVectorRV<V,M>::m_prefix;
  using BaseVectorRV<V,M>::m_imageSet;
  using BaseVectorRV<V,M>::m_pdf;
  using BaseVectorRV<V,M>::m_realizer;
  using BaseVectorRV<V,M>::m_subCdf;
  using BaseVectorRV<V,M>::m_unifiedCdf;
  using BaseVectorRV<V,M>::m_mdf;
};

}  // End namespace QUESO

#endif // UQ_CHEMISTRY_VECTOR_RV_H
