#ifndef UQ_CHEMISTRY_REALIZER_H
#define UQ_CHEMISTRY_REALIZER_H

#include <queso/Defines.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>

namespace QUESO {

/*! 
 * \class ChemistryVectorRealizer
 * \brief A class for handling sampling from a jeffreys probability density
 * distribution.
 *
 * This class handles sampling from a jeffreys probability density distribution.*/

template<class V, class M>
class ChemistryVectorRealizer : public BaseVectorRealizer<V,M> {
public:
  
    //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix and the image set of the vector realizer.  */
  ChemistryVectorRealizer(const char*                  prefix,
                               const VectorSet<V,M>& unifiedImageSet);
  //! Destructor
  ~ChemistryVectorRealizer();
  //@}
  
  //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a jeffreys distribution and saves it in \c nextValues. It 
   * internally finds the minimum and the maximum values of the distribution.
   */  
  void realization(V& nextValues) const;
  //@}
  
private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;
};

}  // End namespace QUESO

#endif // UQ_JEFFREYS_REALIZER_H
