//
// Created by jordan on 31/10/24.
//

#ifndef SAILS_SAILS_MORPH_H
#define SAILS_SAILS_MORPH_H

#include "gemmi/model.hpp" // for Structure

#include "../include/sails-glycan.h"

namespace Sails {

    class Morph {
    public:
        Morph(gemmi::Structure* structure): m_structure(structure) {}

        void transform(Sails::Glycan* glycan, const std::string& WURCS);

    private:
        gemmi::Structure* m_structure;

    };


}
#endif //SAILS_SAILS_MORPH_H
