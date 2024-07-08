//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include "../cpp/sails.cpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(sails_module, m) {
    m.def("run", &run);
}