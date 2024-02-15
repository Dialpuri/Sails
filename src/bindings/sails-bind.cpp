//
// Created by Jordan Dialpuri on 15/02/2024.
//

#include "sails-bind.h"
#include <nanobind/nanobind.h>

int add(int a, int b) { return a + b; }

NB_MODULE(pysails, m) {
m.def("add", &add);
}