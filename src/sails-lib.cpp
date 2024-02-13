//
// Created by Jordan Dialpuri on 12/02/2024.
//

#include "sails-lib.h"

/*
 * SAILS INPUT
 */

/**
 * @brief Load MTZ file and process reflection data
 *
 * Loads the specified MTZ file and processes the reflection data based on the given input columns. It initializes the work reflection data, computes the necessary values, and performs a Fast Fourier Transform (FFT) on the processed data.
 *
 * @param input_mtz The path to the input MTZ file.
 * @param ipcol_fo The column label for the observed structure factors.
 * @param ipcol_hl The column label for the Hendrickson Lattman coefficients.
 * @param ipcol_pw The column label for the phase weight.
 * @param ipcol_fc The column label for the calculated structure factors.
 * @param ipcol_fr The column label for the reflection flags.
 * @param res_in The desired resolution for the processed data.
 *
 * @note This method assumes that the necessary external libraries are available and properly configured.
 * @note The specified input column labels should represent valid columns in the MTZ file.
 *
 * @see load_map()
 * @see load_pdb()
 */
void SailsInput::load_mtz(const std::string& input_mtz, const std::string& ipcol_fo, const std::string& ipcol_hl, const std::string& ipcol_pw, const std::string& ipcol_fc, const std::string& ipcol_fr, float res_in) {

    std::ifstream file(input_mtz);
    if (!file.good()) {

        return;
    }
    using clipper::data32::Compute_fphi_from_fsigf_phifom;
    clipper::Resolution resol;
    mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtzfile.open_read( input_mtz );
    double res = clipper::Util::max( mtzfile.resolution().limit(), clipper::ftype(res_in));

    mtzfile.close_read();
    resol = clipper::Resolution( res );

    // Get work reflection data
    clipper::HKL_info hkls;
    mtzfile.open_read( input_mtz );
    hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );

    clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
    clipper::HKL_data<clipper::data32::ABCD>    wrk_hl( hkls );
    clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw( hkls );
    clipper::HKL_data<clipper::data32::F_phi>   fphi  ( hkls );
    clipper::HKL_data<clipper::data32::Flag>    flag  ( hkls );
    if ( ipcol_fo != "NONE" ) mtzfile.import_hkl_data( wrk_f ,ipcol_fo );
    if ( ipcol_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_hl );
    if ( ipcol_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_pw );
    if ( ipcol_fc != "NONE" ) mtzfile.import_hkl_data( fphi,  ipcol_fc );
    if ( ipcol_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_fr );
    mtzfile.close_read();

    clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;

    for ( clipper::HKL_data_base::HKL_reference_index ih = hkls.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_f1[ih] = clipper::data32::F_sigF();

    clipper::Spacegroup cspg = hkls.spacegroup();
    clipper::Cell cxtl = hkls.cell();
    clipper::Grid_sampling grid = clipper::Grid_sampling( cspg, cxtl, hkls.resolution() );

    xwrk = clipper::Xmap<float>( cspg, cxtl, grid );

    if ( ipcol_hl == "NONE" )
        wrk_hl.compute( wrk_pw, clipper::data32::Compute_abcd_from_phifom() );

    if ( ipcol_pw == "NONE" )
        wrk_pw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );

    if ( ipcol_fc == "NONE" )
        fphi.compute( wrk_f1, wrk_pw, Compute_fphi_from_fsigf_phifom() );

    xwrk.fft_from( fphi );
}

/**
 * \brief Load a CCP4 map file into the SailsInput object.
 *
 * This method opens the CCP4 map file specified by the input_prediction parameter,
 * imports the data to the xpred member variable of type clipper::Xmap<float>,
 * and then closes the file.
 *
 * \param input_prediction The path to the CCP4 map file to load.
 *
 * \return None.
 */
void SailsInput::load_map(const std::string& input_prediction) {
    std::ifstream file(input_prediction);
    if (!file.good()) {
        std::cout << "No suitable map file found, not loading predicted map\n";
        return;
    }

    clipper::CCP4MAPfile map_file;
    map_file.open_read(input_prediction);
    map_file.import_xmap( xpred );
    map_file.close_read();
}

/**
 * @brief Load a PDB file and import it into a MiniMol object
 *
 * @param input_pdb The path to the PDB file to be loaded
 */
void SailsInput::load_pdb(const std::string& input_pdb) {

    std::ifstream file(input_pdb);
    if (!file.good()) {
        std::cout << "No suitable pdb file found, not loading\n";
        return;
    }

    clipper::MMDBfile mfile;
    mfile.read_file(input_pdb);
    mfile.import_minimol(mol);
}

/**
 * @brief Constructs a SailsInput object.
 *
 * This method initializes a SailsInput object by loading and processing input data.
 *
 * @param input_pdb The path to the PDB file.
 * @param input_mtz The path to the MTZ file.
 * @param ipcol_fo The column name for the observed structure factors in the MTZ file.
 * @param ipcol_hl The column name for the Hendrickson Lattman coefficients.
 * @param ipcol_pw The column name for the phase weight in the MTZ file.
 * @param ipcol_fc The column name for the calculated structure factors in the MTZ file.
 * @param ipcol_fr The column name for the reflection flags.
 * @param res_in The resolution to be processed at.
 * @param input_prediction The path to the predicted map file.
 */
SailsInput::SailsInput( const std::string& input_pdb,
                        const std::string& input_mtz,
                        const std::string& ipcol_fo,
                        const std::string& ipcol_hl,
                        const std::string& ipcol_pw,
                        const std::string& ipcol_fc,
                        const std::string& ipcol_fr,
                        float res_in,
                        const std::string& input_prediction) {

    // Create work xmap from MTZ and column names
    load_mtz(input_mtz, ipcol_fo, ipcol_hl, ipcol_pw, ipcol_fc, ipcol_fr, res_in);

    // Get predicted xmap
    load_map(input_prediction);

    // Get PDB model
    load_pdb(input_pdb);
}
/* END */


/*
 * SAILS MONOMERS
 */

/**
 * @brief Load a monomer from a PDB file
 *
 * Loads a monomer from the specified PDB file. It reads the file, imports it into an MMDBfile object, and returns the first monomer from the MiniMol object.
 *
 * @param code The code for the monomer to load.
 *
 * @return The loaded monomer.
 *
 * @note This method assumes that the necessary external libraries are available and properly configured.
 * @note If the specified monomer file cannot be found, the method will output an error message and exit the program.
 *
 * @see load_monomer() in sails-lib.cpp
 * @see MMonomer class in minimol.h
 * @see MMDBfile class in minimol_io.h
 */
clipper::MMonomer SailsMonomers::load_momomer(const std::string& code) {

    std::string path = "data/monomers/" + code + ".pdb";
    std::ifstream file(path);
    if (!file.good()) {
        std::cout << "No suitable monomer file found at " << path << ", aborting\n";
        exit(-1);
    }

    clipper::MMDBfile mfile;
    clipper::MiniMol mol;
    mfile.read_file(path);
    mfile.import_minimol(mol);

    return mol[0][0];
}

/* END */
