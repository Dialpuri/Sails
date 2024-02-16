//
// Created by Jordan Dialpuri on 14/02/2024.
//
#pragma once

#ifndef SAILS_LIB_H
#define SAILS_LIB_H

#include <clipper/clipper.h>
#include "sails-model.h"

namespace Sails {

    class TargetFunctionBase {
    public:
        TargetFunctionBase() {}

        virtual ~TargetFunctionBase() {}

        virtual int num_params() const = 0;

        virtual double operator()(const std::vector<double> &args) const = 0;
    };

    class SimplexOptimiser {
    public:
        enum TYPE {
            NORMAL, GRADIENT
        };

        SimplexOptimiser(double tolerance = 1e-3, int max_cycles = 50, TYPE type = NORMAL, bool debug = false) :
                m_tolerance(tolerance), m_max_cycles(max_cycles), m_type(type), m_debug_mode(debug) {};

        std::vector<double> operator()(const TargetFunctionBase &target_function, const std::vector<std::vector<double>> &args) const;

        void debug() { m_debug_mode = true; }

    private:
        double m_tolerance;
        double m_max_cycles;
        TYPE m_type;
        bool m_debug_mode;
    };

    class TargetFunctionSugar : public TargetFunctionBase {
    public:
        typedef float (*Scorer)(Sails::LinkageSet*, clipper::Xmap<float> &);

        TargetFunctionSugar() = default;

        TargetFunctionSugar(clipper::Xmap<float> &xmap,
                            Sails::LinkageSet* linkage,
                            Scorer scorer,
                            float step
        ) :
                m_xmap(&xmap),
                m_linkage(linkage),
                m_scorer(scorer),
                m_step(step) {};

        int num_params() const override { return 3; }

        double operator()(const std::vector<double> &args) const override;

        LinkageSet * refine();

    private:
        clipper::Xmap<float> *m_xmap{};
        Sails::LinkageSet* m_linkage;
        Scorer m_scorer{};
        float m_step = clipper::Util::d2rad(1);
    };
}
#endif 