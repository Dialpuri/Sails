#include "sails-refine.h"

std::vector<double> Sails::SimplexOptimiser::operator()(const TargetFunctionBase& target_function, const std::vector<std::vector<double>>& args) const { 
    if (m_debug_mode) std::cout << "() operator called" << std::endl;

    enum STEP { UNKN, EXTN, NRML, CTRN, CTRX };

    // check input
    int size = target_function.num_params();
    if ( args.size() != size+1 )
        clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );
    for ( int i = 0; i < args.size(); i++ )
        if ( args[i].size() != size )
            clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );

    // make arrays
    std::vector<std::vector<double> > params( size + 1 );
    std::vector<double> fn( size + 1 );
    // calc target fn
    for ( int i = 0; i < args.size(); i++ ) {
        params[i] = args[i];
        fn[i] = target_function( args[i] );
    }

    // simplex loop
    int worst_index, best_index;
    worst_index = best_index = 0;
    double f0(0.0), f1(0.0), f2(0.0);
    for ( int cyc = 0; cyc < m_max_cycles; cyc++ ) {
        if (m_debug_mode)
            std::cout << cyc << "/" << m_max_cycles << std::endl;

        if ( m_debug_mode )  // DEBUG OUTPUT
            for ( int i = 0; i < params.size(); i++ ) {
                std::cout << "Index " << i << "  score of i " << clipper::String( fn[i] ) << "  ";
                for ( int j = 0; j < size; j++ )
                    std::cout << "params[" << i << "][" << j << "] " << clipper::String( params[i][j] ) << "\t";
                std::cout << "\n";
            }
        // find worst point
        worst_index = best_index = 0;
        for ( int i = 0; i < params.size(); i++ ) {
            if ( fn[i] > fn[worst_index] ) worst_index = i;
            if ( fn[i] < fn[best_index] ) best_index = i;
        }

        if (m_debug_mode)
            std::cout << "worst " << worst_index << "\tbest " << best_index << std::endl;

        // termination condition
        if (fn[worst_index] - fn[best_index] < m_tolerance )
        {
            if (m_debug_mode)
                std::cout << "Termination hit " << std::endl;
            break;
        }
        // find centroid of the rest
        std::vector<double> centr( size, 0.0 ),
                shift( size ), t0( size ), t1( size ), t2( size );
        double sumweight = 0.0;

        for ( int i = 0; i < params.size(); i++ )
            if (i != worst_index ) {
                double weight = 1.0;
                if ( m_type == GRADIENT ) {
                    double r2 = 0.0;
                    for ( int j = 0; j < size; j++ )
                        r2 += clipper::Util::sqr( params[i][j] - params[worst_index][j] );
                    weight = (fn[worst_index] - fn[i] ) / sqrt(r2 );
                }
                for ( int j = 0; j < size; j++ )
                    centr[j] += weight * params[i][j];
                sumweight += weight;
            }

        for ( int j = 0; j < size; j++ ) centr[j] = centr[j] / sumweight;
        for ( int j = 0; j < size; j++ ) shift[j] = centr[j] - params[worst_index][j];
        // calculate first trial point
        for ( int j = 0; j < size; j++ ) {
//            std::cout << centr[j] << "\t" << shift[j] << std::endl;
            t0[j] = centr[j] - 0.5 * shift[j];
            t1[j] = centr[j] + 1.0 * shift[j];
            t2[j] = centr[j] + 2.0 * shift[j];
        }

        f1 = target_function( t1 );
        // simplex conditions
        STEP step = UNKN;

        if ( !clipper::Util::is_nan( f1 ) ) { // new point is valid
            if ( f1 < fn[worst_index] ) {    // new point is better than worst
                if ( f1 < fn[best_index] ) {  // new point is better than best
                    f2 = target_function( t2 );
                    if ( !clipper::Util::is_nan( f2 ) ) {
                        if ( f2 < f1 ) {  // extended step is best
                            step = EXTN;
                        } else {  // normal step better than extended step
                            step = NRML;
                        }
                    } else {  // extended step is invalid
                        step = NRML;
                    }
                } else {  // normal step is neither best nor worst (default)
                    step = NRML;
                }
            }  // else normal step is worse
        }  // else normal step is invalid
        if ( step == UNKN ) {    // normal step is worse or invalid
            f0 = target_function( t0 );  // contraction step (always valid)
            if ( f0 < fn[worst_index] ) {   // contraction step is better than worst
                step = CTRN;
            } else {               // contraction step is worse
                step = CTRX;
            }
        }
        if      ( step == EXTN ) { params[worst_index] = t2; fn[worst_index] = f2; }  // extension
        else if ( step == NRML ) { params[worst_index] = t1; fn[worst_index] = f1; }  // normal
        else if ( step == CTRN ) { params[worst_index] = t0; fn[worst_index] = f0; }  // contraction
        else {  // otherwise contract all towards best (slow)
            for ( int i = 0; i < params.size(); i++ )
                if (i != best_index ) {
                    for ( int j = 0; j < size; j++ )
                        params[i][j] = 0.5 * ( params[i][j] + params[best_index][j] );
                    fn[i] = target_function( params[i] );
                }
        }
        if ( m_debug_mode ) {  // DEBUG OUTPUT
            if      ( step == EXTN ) std::cout << "Extn-step\n";
            else if ( step == NRML ) std::cout << "Nrml-step\n";
            else if ( step == CTRN ) std::cout << "Ctrn-step\n";
            else                     std::cout << "Ctrx-step\n";
        }
    }
    return params[best_index];
}

double Sails::TargetFunctionSugar::operator()(const std::vector<double>& args) const { 
    clipper::Euler_ccp4 euler = {args[0], args[1], args[2]};
    clipper::Rotation rotation(euler);
    
    clipper::RTop_orth rtop = {rotation.matrix(), m_translation};
    clipper::MMonomer test_fragment = m_sugar;
    test_fragment.transform(rtop);

    float score = m_scorer(test_fragment, *m_xmap);
    return -score; 
}

clipper::RTop_orth Sails::TargetFunctionSugar::refine(Sails::LinkageSet& linkage) { 
    std::vector<double> args = {0,0,0};
    std::vector<std::vector<double>> args_init = {
        args, 
        {args[0]+m_step, args[1], args[2]},
        {args[0], args[1]+m_step, args[2]},
        {args[0], args[1], args[2]+m_step}
    }; 

    double tol = 5e-4 * (*this)(args);
    SimplexOptimiser optimiser = SimplexOptimiser(tol, 20, SimplexOptimiser::NORMAL, false);
    std::vector<double> refined_arguments = optimiser(*this, args_init);

    clipper::Euler_ccp4 refined_euler = clipper::Euler_ccp4(refined_arguments[0], refined_arguments[1], refined_arguments[2]);
    clipper::Rotation refined_rotation = clipper::Rotation(refined_euler);

    return {refined_rotation.matrix(), m_translation};
}