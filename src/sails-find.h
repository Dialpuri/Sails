//
// Created by Jordan Dialpuri on 12/02/2024.
//

#ifndef SAILS_FIND_H
#define SAILS_FIND_H

#include <clipper/clipper-minimol.h>
#include "sails-lib.h"
#include "sails-model.h"
#include "sails-util.h"


class SailsFind {
public:
    SailsFind(clipper::MiniMol& work_model);

    clipper::MiniMol find();

private:
    clipper::MiniMol mol;
};
#endif //SAILS_FIND_H
