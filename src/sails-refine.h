//
// Created by Jordan Dialpuri on 14/02/2024.
//

#ifndef SAILS_LIB_H
#define SAILS_LIB_H

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>
#include "sails-model.h"

namespace Sails { 

    class TargetFunctionBase { 
    public:
        TargetFunctionBase() {}
        virtual ~TargetFunctionBase() {}
        virtual int num_params() const = 0; 
        virtual double operator()(const std::vector<double>& args) const = 0; 
    };

    class SimplexOptimiser { 
    public: 
        enum TYPE { 
            NORMAL, GRADIENT
        };
    
        SimplexOptimiser(double tolerance = 1e-3, int max_cycles = 50, TYPE type = NORMAL, bool debug = false): 
        m_tolerance(tolerance), m_max_cycles(max_cycles), m_type(type), m_debug_mode(debug) {};

        std::vector<double> operator()(const TargetFunctionBase& target_function, const std::vector<std::vector<double>>& args) const;

        void debug() {m_debug_mode = true;}

    private:
        double m_tolerance;
        double m_max_cycles;
        TYPE m_type; 
        bool m_debug_mode;
    };

    class TargetFunctionSugar: TargetFunctionBase { 
    public: 
        typedef float (*Scorer)(clipper::MMonomer&, clipper::Xmap<float>&);
        TargetFunctionSugar() = default;
        TargetFunctionSugar(clipper::Xmap<float>& xmap,
                            clipper::Vec3<> translation,
                            clipper::MMonomer sugar, 
                            Scorer scorer, 
                            float step                           
                            ): 
                            m_xmap(&xmap), 
                            m_step(step),
                            m_translation(translation),
                            m_scorer(scorer) {};

        inline int num_params() const override {return 3; }

        double operator()(const std::vector<double>& args) const; 

        clipper::RTop_orth refine(Sails::LinkageSet& linkage); 

    private:
        clipper::Xmap<float> *m_xmap{};
        clipper::Vec3<> m_translation; 
        clipper::MMonomer m_sugar; 
        Scorer m_scorer;
        float m_step = clipper::Util::d2rad(1);

    };


    class Refine { 
        Refine(clipper::MiniMol& mol);
    };

}

#endif 