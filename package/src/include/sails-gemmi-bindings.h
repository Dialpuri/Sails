//
// Created by Jordan Dialpuri on 22/07/2024.
//

#ifndef SAILS_OUTPUT_H
#define SAILS_OUTPUT_H

#include <string>
#include <iostream>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/mtz.hpp>
#include <gemmi/unitcell.hpp>

namespace Sails {
    struct NumberPair {
        NumberPair() = default;
        NumberPair(double value1, double value2)
            : value1(value1),
              value2(value2) {
        }

        double value1 = std::numeric_limits<double>::quiet_NaN();
        double value2 = std::numeric_limits<double>::quiet_NaN();
    };

    struct HKL {
        HKL(int h, int k, int l)
            : h(h),
              k(k),
              l(l) {
        }

        int h, k, l;
    };

    struct Reflection {
        Reflection(const HKL &hkl, const NumberPair &f_sigf)
                  : hkl(hkl),
                    f_sigf(f_sigf) {
        }

        Reflection(const HKL &hkl, const NumberPair &f_sigf, const NumberPair &fwt_phwt,
                   const NumberPair &delfwt_phdelwt)
                   : hkl(hkl),
                     f_sigf(f_sigf),
                     fwt_phwt(fwt_phwt),
                     delfwt_phdelwt(delfwt_phdelwt) {
        }

        HKL hkl;
        NumberPair f_sigf;
        NumberPair fwt_phwt;
        NumberPair delfwt_phdelwt;
    };

    struct Cell {
        explicit Cell(const gemmi::UnitCell& cell):
                    a(cell.a), b(cell.b), c(cell.c), alpha(cell.alpha), beta(cell.beta), gamma(cell.gamma) {}

        Cell(double a, double b, double c, double alpha, double beta, double gamma):
            a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {}

        double a, b, c, alpha, beta, gamma;

        [[nodiscard]] gemmi::UnitCell to_gemmi_cell() const {return {a, b, c, alpha, beta, gamma};}
    };

    struct MTZ {
        MTZ(std::vector<Reflection>& reflections, Cell& cell, std::string& spacegroup):
            reflections(reflections), cell(cell), spacegroup(spacegroup) {}

        std::vector<Reflection> reflections;
        Cell cell;
        std::string spacegroup;
    };



    gemmi::Mtz form_gemmi_mtz(MTZ& mtz);

    MTZ form_sails_mtz(gemmi::Mtz& mtz);


    struct Output {
        Output(gemmi::Structure& structure, MTZ& mtz): structure(structure), mtz(mtz) {};

        gemmi::Structure structure ;
        MTZ mtz;
    };

}

#endif //SAILS_OUTPUT_H
