//
// Created by Jordan Dialpuri on 12/02/2024.
//

#ifndef SAILS_LIB_H
#define SAILS_LIB_H

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <fstream>

class SailsInput {
public:

    explicit SailsInput(const std::string& input_pdb,
                        const std::string& input_mtz,
                        const std::string& ipcol_fo,
                        const std::string& ipcol_hl,
                        const std::string& ipcol_pw,
                        const std::string& ipcol_fc,
                        const std::string& ipcol_fr,
                        float res_in,
                        const std::string& input_prediction);

    clipper::Xmap<float> get_work_xmap() const noexcept {
        return xwrk;
    }

    clipper::Xmap<float> get_predicted_map() const noexcept {
        return xpred;
    }

    clipper::MiniMol get_minimol() const noexcept {
        return mol;
    }
private:
    void load_mtz(const std::string&input_mtz,
                const std::string& ipcol_fo,
                const std::string& ipcol_hl,
                 const std::string& ipcol_pw,
                 const std::string& ipcol_fc,
                 const std::string& ipcol_fr,
                 float res_in);

    void load_map(const std::string& input_prediction);
    void load_pdb(const std::string& input_pdb);

private:
    clipper::Xmap<float> xwrk;
    clipper::Xmap<float> xpred;
    clipper::CCP4MTZfile mtzfile;
    clipper::MiniMol mol;

};

class SailsMonomers {
public:
    SailsMonomers();

    static clipper::MMonomer load_momomer(const std::string& code);
};
#endif //SAILS_LIB_H
