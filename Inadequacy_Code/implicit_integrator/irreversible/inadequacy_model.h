#ifndef INADEQUACY_MODEL_H
#define INADEQUACY_MODEL_H

#include <vector>

class inadequacy_model
{
  public:
    inadequacy_model (unsigned int n_inad_from_user, 
                      unsigned int n_reactions_inad_from_user);
    unsigned int n_inad;            // number of virtual species
    unsigned int n_reactions_inad;  // number of virtual species
};

#endif
