
#include "sails-find.h"
#include "sails-lib.h"

void run_sails( const std::string& input_pdb,
                const std::string& input_mtz,
                const std::string& ipcol_fo,
                const std::string& ipcol_hl,
                const std::string& ipcol_pw,
                const std::string& ipcol_fc,
                const std::string& ipcol_fr,
                float res_in,
                const std::string& input_prediction,
                const std::string& output_pdb) {
    SailsInput sails_input = SailsInput(input_pdb, input_mtz, ipcol_fo, ipcol_hl, ipcol_pw, ipcol_fc, ipcol_fr, res_in, input_prediction);
    clipper::Xmap<float> work_map = sails_input.get_work_xmap();
    clipper::Xmap<float> pred_map = sails_input.get_predicted_map();
    clipper::MiniMol work_model = sails_input.get_minimol();

    SailsFind sails_find = SailsFind(work_model, work_map, pred_map);
    work_model = sails_find.find();

    clipper::MMDBfile mf;
    mf.export_minimol(work_model);
    mf.write_file(output_pdb, clipper::MMDBfile::TYPE::PDB);
}

int main(int argc, char **argv) {
    CCP4Program prog("sails", "0.0.1", "$Date: 2024/2/12");
    prog.set_termination_message("Failed");

    std::cout << std::endl << "Copyright 2024 Jordan Dialpuri and University of York." << std::endl << std::endl;
    prog.summary_beg();
    prog.summary_end();

    clipper::String input_pdb = "";
    clipper::String input_mtz = "";
    clipper::String input_prediction = "";

    clipper::String ipcol_fo =  "NONE";
    clipper::String ipcol_hl =  "NONE";
    clipper::String ipcol_pw =  "NONE";
    clipper::String ipcol_fc =  "NONE";
    clipper::String ipcol_fr =  "NONE";
    double res_in = 1.0;

    CCP4CommandInput args(argc, argv, true);
    int arg = 0;
    while (++arg < args.size()) {
        if (args[arg] == "-pdbin") {
            if (++arg < args.size())
                input_pdb = args[arg];
        }
        else if (args[arg] == "-mtzin") {
            if (++arg < args.size())
                input_mtz = args[arg];
        }
        else if ( args[arg] == "-colin-fo" )
        {
            if ( ++arg < args.size() )
                ipcol_fo = args[arg];
        }
        else if ( args[arg] == "-colin-hl" )
        {
            if ( ++arg < args.size() )
                ipcol_hl = args[arg];
        }
        else if ( args[arg] == "-colin-phifom" )
        {
            if ( ++arg < args.size() )
                ipcol_pw = args[arg];
        }
        else if ( args[arg] == "-colin-fc" )
        {
            if ( ++arg < args.size() )
                ipcol_fc = args[arg];
        }
        else if ( args[arg] == "-colin-free" )
        {
            if ( ++arg < args.size() )
                ipcol_fr = args[arg];
        }
        else if ( args[arg] == "-resolution" )
        {
            if ( ++arg < args.size() )
                res_in = clipper::String(args[arg]).f();
        }
        else if ( args[arg] == "-predin" )
        {
            if ( ++arg < args.size() )
                input_prediction = args[arg];
        }
        else {
            std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    SailsInput sails_input = SailsInput(input_pdb, input_mtz, ipcol_fo, ipcol_hl, ipcol_pw, ipcol_fc, ipcol_fr, res_in, input_prediction);
    clipper::Xmap<float> work_map = sails_input.get_work_xmap();
    clipper::Xmap<float> pred_map = sails_input.get_predicted_map();
    clipper::MiniMol work_model = sails_input.get_minimol();

    SailsFind sails_find = SailsFind(work_model, work_map, pred_map);
    work_model = sails_find.find();


    prog.set_termination_message( "Normal termination" );
    return 0;
}

