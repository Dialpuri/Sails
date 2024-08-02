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

    /**
     * @brief Represents a pair of numbers which default to NaN
     *
     * Used to represent a set of floating point numbers, such as F and SIGF reflection values
     */
    struct NumberPair {
        NumberPair() = default;
        NumberPair(double value1, double value2)
            : value1(value1),
              value2(value2) {
        }

        double value1 = std::numeric_limits<double>::quiet_NaN();
        double value2 = std::numeric_limits<double>::quiet_NaN();
    };

    /**
     * @brief Represents the Miller indices (h, k, l) of a reflection
     *
     * This struct is used to represent the Miller indices (h, k, l) of a reflection in crystallography.
     * It provides a way to store and access these indices as individual integer values.
     */
    struct HKL {
        HKL(int h, int k, int l)
            : h(h),
              k(k),
              l(l) {
        }

        int h, k, l;
    };

    /**
     * @brief Reflection represents a reflection in a crystallographic dataset
     *
     * The Reflection class is used to store various properties of a reflection, including its HKL values and associated
     * numerical pairs. It provides constructors to initialize different combinations of properties.
     */
    struct Reflection {
        Reflection(const HKL &hkl, const NumberPair &f_sigf)
                  : hkl(hkl),
                    f_sigf(f_sigf) {
        }

        Reflection(const HKL &hkl, const NumberPair &f_sigf, const NumberPair &fwt_phwt)
                  : hkl(hkl),
                    f_sigf(f_sigf),
                    fwt_phwt(fwt_phwt){
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

    /**
     * @brief Represents a unit cell in three-dimensional space
     *
     * This struct provides a representation of a unit cell, which is a parallelepiped in three-dimensional space.
     * It stores the parameters that define the unit cell: a, b, c (lengths of the edges) and alpha, beta, gamma (angles between the edges).
     * It also provides a method to convert the Cell object to a gemmi::UnitCell object.
     */
    struct Cell {
        explicit Cell(const gemmi::UnitCell& cell):
                    a(cell.a), b(cell.b), c(cell.c), alpha(cell.alpha), beta(cell.beta), gamma(cell.gamma) {}

        Cell(double a, double b, double c, double alpha, double beta, double gamma):
            a(a), b(b), c(c), alpha(alpha), beta(beta), gamma(gamma) {}

        double a, b, c, alpha, beta, gamma;

        [[nodiscard]] gemmi::UnitCell to_gemmi_cell() const {return {a, b, c, alpha, beta, gamma};}
    };

    /**
     * @brief Represents a crystallographic MTZ file within Sails
     *
     * This class is used to store information about a crystallographic MTZ file,
     * which is a common file format used in X-ray crystallography to store experimental data. It contains
     * a vector of Reflection objects representing the observed reflections, a Cell object representing the
     * unit cell parameters, and a string representing the space group.
     */
    struct MTZ {
        MTZ(std::vector<Reflection>& reflections, Cell& cell, std::string& spacegroup):
            reflections(reflections), cell(cell), spacegroup(spacegroup) {}

        std::vector<Reflection> reflections;
        Cell cell;
        std::string spacegroup;
    };


    /**
     * @brief Converts a Sails::MTZ to a gemmi::Mtz object
     *
     * This function takes an Sails::MTZ and converts it to a gemmi::Mtz object.
     * It extracts the necessary information from the Sails::MTZ, such as
     * unit cell, space group, reflection data, and creates a new gemmi::Mtz object
     * with the converted data.
     *
     * @param mtz The input MTZ file to be converted
     * @return A gemmi::Mtz object representing the converted MTZ*/
    gemmi::Mtz form_gemmi_mtz(MTZ& mtz);

    /**
     * @brief Forms a Sails::MTZ from a given gemmi::Mtz object
     *
     * This function takes a gemmi::Mtz object and extracts the necessary columns to form a Sails MTZ representation.
     * It checks the column sizes and throws an exception if they are unequal. It then iterates over the columns
     * and forms a Sails Reflection object for each reflection. Finally, it creates a Sails MTZ object and returns it.
     *
     * @param mtz The gemmi::Mtz object containing the necessary columns
     * @param f_label The F column name, defaults to "FP"
     * @param sigf_label The SIGF column name , defaults to "SIGF"
     * @return The Sails MTZ object formed from the gemmi::Mtz object
     */
    MTZ form_sails_mtz(gemmi::Mtz& mtz, const std::string &f_label="F", const std::string &sigf_label = "SIGF");


    /**
     * @brief Represents the output of a Sails process
     *
     * This class is used to hold the output of a process, which includes a gemmi::Structure object, an MTZ object and
     * log string.
     */
    struct Output {
        Output(gemmi::Structure& structure, MTZ& mtz, std::string log): structure(structure), mtz(mtz), log(std::move(log)) {};

        gemmi::Structure structure ;
        MTZ mtz;
        std::string log;
    };

}

#endif //SAILS_OUTPUT_H
