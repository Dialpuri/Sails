//
// Created by jordan on 05/07/24.
//

#include "../include/sails-json.h"


Sails::JSONLoader::JSONLoader(const std::string &file_path) {
    this->init(file_path);
}

void Sails::JSONLoader::init(const std::string &file_path) {
    if (!Utils::file_exists(file_path)) { throw std::runtime_error("Data file could not be found"); }
    m_path = file_path;
}


Sails::ResidueDatabase Sails::JSONLoader::load_residue_database() {

    simdjson::ondemand::parser parser;
    auto json = simdjson::padded_string::load(m_path);
    m_doc = parser.iterate(json);

    // load residues
    const char *residue_key = "residues";
    const char *name_key = "name";
    const char *donor_set_key = "donorSets";
    const char *acceptor_set_key = "acceptorSets";
    const char *snfg_shape_key = "snfgShape";
    const char *snfg_colour_key = "snfgColour";
    const char *preferred_depth_key = "preferredDepths";
    const char *anomer_key = "anomer";
    const char *wurcs_code_key = "wurcsCode";
    const char *special_key = "special";


    ResidueDatabase database;
    auto residues = m_doc[residue_key];
    for (auto value: residues) {
        std::string name = std::string(value[name_key].get_string().value());

        std::vector<AtomSet> donor_sets = extract_atom_set(value, donor_set_key);
        std::vector<AtomSet> acceptors_sets = extract_atom_set(value, acceptor_set_key);

        std::string snfg_shape = std::string(value[snfg_shape_key].get_string().value());
        std::string snfg_colour = std::string(value[snfg_colour_key].get_string().value());

        std::vector<int> preferred_depths;
        for (auto a: value[preferred_depth_key]) {
            preferred_depths.emplace_back(static_cast<int>(a.get_int64()));
        }

        std::string anomer = std::string(value[anomer_key].get_string().value());
        std::string wurcs_code = std::string(value[wurcs_code_key].get_string().value());

        bool special = value[special_key].get_bool();
        ResidueData data = {acceptors_sets, donor_sets, snfg_shape, snfg_colour, preferred_depths, anomer, wurcs_code, special};
        database.insert({name, data});
    }

    return database;
}

Sails::LinkageDatabase Sails::JSONLoader::load_linkage_database() {
    simdjson::ondemand::parser parser;
    auto json = simdjson::padded_string::load(m_path);
    m_doc = parser.iterate(json);

    // load linkages
    const char *linkage_key = "linkages";
    const char *donor_residue_key = "donorResidue";
    const char *acceptor_residue_key = "acceptorResidue";
    const char *donor_number_key = "donorNumber";
    const char *acceptor_number_key = "acceptorNumber";
    const char *length_key = "length";
    const char *angles_key = "angles";
    const char *torsions_key = "torsions";
    const char *cluster_key = "clusters";

    LinkageDatabase database;

    auto linkages = m_doc[linkage_key];

    for (auto value: linkages) {
        auto donor_residue = std::string(value[donor_residue_key].get_string().value());
        auto acceptor_residue = std::string(value[acceptor_residue_key].get_string().value());
        auto donor_number = int(value[donor_number_key].get_number().value().get_int64());
        auto acceptor_number = int(value[acceptor_number_key].get_number().value().get_int64());

        auto length = value[length_key].get_number().value().as_double();

        Clusters clusters = {};
        for (auto cluster: value[cluster_key]) {
            auto angles = cluster[angles_key];
            AngleSet angle_set = extract_angles(angles);

            auto torsions = cluster[torsions_key];
            TorsionSet torsion_set = extract_torsions(torsions);

            clusters.emplace_back(angle_set, torsion_set);
        }

        LinkageData data = {donor_residue, acceptor_residue, donor_number, acceptor_number, length, clusters};
        database[donor_residue].emplace_back(data);
    }

    return database;
}

Sails::AngleSet Sails::JSONLoader::extract_angles(simdjson::simdjson_result<simdjson::ondemand::value> &value) {
    const char *alpha_mean_key = "alphaMean";
    const char *alpha_stddevkey = "alphaStdDev";
    const char *beta_mean_key = "betaMean";
    const char *beta_stddevkey = "betaStdDev";
    const char *gamma_mean_key = "gammaMean";
    const char *gamma_stddev_key = "gammaStdDev";

    auto alpha = Angle(value[alpha_mean_key].get_number().value().as_double(),
                       value[alpha_stddevkey].get_number().value().as_double());
    auto beta = Angle(value[beta_mean_key].get_number().value().as_double(),
                      value[beta_stddevkey].get_number().value().as_double());
    auto gamma = Angle(value[gamma_mean_key].get_number().value().as_double(),
                       value[gamma_stddev_key].get_number().value().as_double());

    return {alpha, beta, gamma};
}


Sails::TorsionSet Sails::JSONLoader::extract_torsions(simdjson::simdjson_result<simdjson::ondemand::value> &value) {
    const char *phi_mean_key = "phiMean";
    const char *phi_stddev_key = "phiStdDev";
    const char *psi_mean_key = "psiMean";
    const char *psi_stddev_key = "psiStdDev";
    const char *omega_mean_key = "omegaMean";
    const char *omega_stddev_key = "omegaStdDev";

    auto phi = Angle(value[phi_mean_key].get_number().value().as_double(),
                     value[phi_stddev_key].get_number().value().as_double());
    auto psi = Angle(value[psi_mean_key].get_number().value().as_double(),
                     value[psi_stddev_key].get_number().value().as_double());
    auto omega = Angle(value[omega_mean_key].get_number().value().as_double(),
                       value[omega_stddev_key].get_number().value().as_double());

    return {psi, phi, omega};
}

void Sails::JSONWriter::write_json_file(TelemetryLog &log, std::ostream &stream) {
    const auto now = std::chrono::system_clock::now();
    const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
    stream << "{\n";
    stream << "\t\"date\": \"" << strtok(ctime(&t_c), "\n") << "\",\n";
    stream << "\t\"cycles\":[\n\t\t";
    for (const auto &[cycle, entries]: log) {
        stream << "{\n";
        stream << "\t\t\t\"cycle\": " << cycle << ",\n";
        stream << "\t\t\t\"entries\": {\n";
        for (int i = 0; i < entries.size(); ++i) {
            stream << "\t\t\t\t\"" << entries[i].residue_id << "\": {\"rscc\": " << entries[i].rscc_score <<
                    ", \"rsr\": " << entries[i].rsr_score <<
                    ", \"dds\": " << entries[i].dds_score << "}";
            if (i < entries.size() - 1) {
                stream << ",";
            }
            stream << "\n";
        }
        stream << "\t\t\t}\n\t\t}";
        if (cycle < log.size()) stream << ",";
    }
    stream << "]\n}";
}


std::vector<Sails::AtomSet>
Sails::JSONLoader::extract_atom_set(simdjson::simdjson_result<simdjson::ondemand::value> &value, const char *key) {
    const char *atom1_key = "atom1";
    const char *atom2_key = "atom2";
    const char *atom3_key = "atom3";
    const char *identifier_key = "identifier";

    simdjson::simdjson_result<simdjson::ondemand::array> sets = value[key].get_array();

    std::vector<Sails::AtomSet> atom_sets;
    for (auto set: sets) {
        std::string atom_1 = std::string(set[atom1_key].get_string().value());
        std::string atom_2 = std::string(set[atom2_key].get_string().value());
        std::string atom_3 = std::string(set[atom3_key].get_string().value());
        int identifier = int(set[identifier_key].get_int64());

        Sails::AtomSet atom_set = {atom_1, atom_2, atom_3, identifier};
        atom_sets.emplace_back(atom_set);
    }

    std::sort(atom_sets.begin(), atom_sets.end(), [](const Sails::AtomSet &a, const Sails::AtomSet &b) {
        return a.identifier < b.identifier;
    });
    return atom_sets;}
