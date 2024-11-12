//
// Created by Jordan Dialpuri on 01/11/2024.
//

#include <emscripten/bind.h>

using namespace emscripten;

int test() {
  return 4;
}

EMSCRIPTEN_BINDINGS(sails_module)
{
    function("test", &test);
}
