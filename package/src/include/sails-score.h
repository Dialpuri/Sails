//
// Created by Jordan Dialpuri on 22/10/2025.
//

#ifndef SAILS_SCORE_H
#define SAILS_SCORE_H
#include "sails-model.h"
#include "density/sails-density.h"

namespace Sails::Score {

    std::map<Glycosite, double> calculate_rsccs(Sails::Density* density, gemmi::Structure* structure, ResidueDatabase &residue_database);

}

#endif //SAILS_SCORE_H
