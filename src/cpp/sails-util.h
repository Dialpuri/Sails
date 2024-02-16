//
// Created by Jordan Dialpuri on 12/02/2024.
//

#ifndef SAILS_UTIL_H
#define SAILS_UTIL_H

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

struct SailsUtil {

    static clipper::MAtom create_atom(const clipper::Coord_orth& coord, std::string id, std::string element = "C") {
        clipper::MAtom m_atom;
        m_atom.set_u_iso ( 0.25 );
        m_atom.set_u_aniso_orth(clipper::U_aniso_orth(clipper::Mat33sym<>::null()));
        m_atom.set_occupancy( 1.0 );
        m_atom.set_id(std::move(id));
        m_atom.set_element( std::move(element));
        m_atom.set_coord_orth(coord);
        return m_atom;
    }
};



#endif //SAILS_UTIL_H
