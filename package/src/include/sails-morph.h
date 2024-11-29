//
// Created by jordan on 31/10/24.
//

#ifndef SAILS_SAILS_MORPH_H
#define SAILS_SAILS_MORPH_H

#include "gemmi/model.hpp" // for Structure
#include "sails-glycan.h"

namespace Sails {
    class Morph {
    public:
        Morph(gemmi::Structure *structure): m_structure(structure) {
        }

        void transform(Glycan &glycan, PseudoGlycan &pseudo_glycan);

        void swap_sugars(Glycan& glycan, PseudoGlycan& pseudo_glycan) const;

        static bool check_graph_isomorphism(Glycan &glycan, PseudoGlycan &pseudo_glycan);

    private:
        gemmi::Structure *m_structure;
    };
}
#endif //SAILS_SAILS_MORPH_H
