//
// Created by Jordan Dialpuri on 16/02/2024.
//
#include <nanobind/nanobind.h>
#include "../cpp/csails.cpp"
#include <nanobind/stl/string.h>

namespace nb = nanobind;

using namespace nb::literals;

NB_MODULE(sails_module, m) {
nb::class_<SailsInput>(m, "Input")
            .def(nb::init< const std::string&, const std::string&, const std::string&, const std::string&,
                    const std::string&, const std::string&, const std::string&, float, const std::string&
                 >());

nb::class_<SailsOutput>(m, "Output")
        .def(nb::init<const std::string&>());

nb::class_<SailsFind>(m, "Find")
        .def(nb::init<SailsInput& >())
        .def(nb::init<SailsInput& , SailsOutput& >());

m.def("run", &run_sails);
}