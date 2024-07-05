//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include "../include/sails-json.h"
#include "../include/sails-model.h"

namespace nb = nanobind;
using namespace nb::literals;

void run() {
    Sails::ResidueDatabase database;
    Sails::JSONLoader loader = {"/home/jordan/dev/Sails/package/data/data.json", database};

}

NB_MODULE(sails_module, m) {
m.def("run", &run);
}