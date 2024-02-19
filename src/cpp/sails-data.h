//
// Created by Jordan Dialpuri on 19/02/2024.
//

#ifndef SAILS_DATA_H
#define SAILS_DATA_H
#include "sails-model.h"

struct Data {
    std::set<std::string> possible_links;
    std::set<std::string> donors;
    std::set<std::string> acceptors;
    Sails::ResidueType type;

    bool can_link(const std::string& residue_name) {
        return std::find(possible_links.begin(), possible_links.end(), residue_name) == possible_links.end();
    }

    [[nodiscard]] virtual std::set<std::string> get_donors() const { return donors;}
    [[nodiscard]] virtual std::set<std::string> get_acceptors() const { return acceptors;}
    [[nodiscard]] virtual std::set<std::string> get_possible_links() const { return possible_links;}

    [[nodiscard]] std::string format() {
        return "Donors: " + std::to_string(get_donors().size()) + " Acceptors: " + std::to_string(get_acceptors().size())
        + " Possible Links " + std::to_string(get_possible_links().size());
    }
};

struct ASNData: Data {
    std::set<std::string> donors = {"ND2"};
    std::set<std::string> acceptors = {};
    std::set<std::string> possible_links = {"NAG", "NGA"};
    Sails::ResidueType type = Sails::ResidueType::ASN;
    [[nodiscard]] std::set<std::string> get_donors() const override { return donors;}
    [[nodiscard]] std::set<std::string> get_acceptors() const override { return acceptors;}
    [[nodiscard]] std::set<std::string> get_possible_links() const override { return possible_links;}
};

struct NAGData: Data {
    std::set<std::string> donors = {"O4", "O6"};
    std::set<std::string> acceptors = {"C1"};
    std::set<std::string> possible_links = {"NAG", "BMA", "MAN", "FUC"};
    Sails::ResidueType type = Sails::ResidueType::NAG;
    [[nodiscard]] std::set<std::string> get_donors() const override { return donors;}
    [[nodiscard]] std::set<std::string> get_acceptors() const override { return acceptors;}
    [[nodiscard]] std::set<std::string> get_possible_links() const override { return possible_links;}
};


struct BMAData: Data {
    std::set<std::string> donors = {"O3", "O6"};
    std::set<std::string> acceptors = {"C1"};
    std::set<std::string> possible_links = {"NAG", "BMA", "MAN", };
    Sails::ResidueType type = Sails::ResidueType::NAG;
    [[nodiscard]] std::set<std::string> get_donors() const override { return donors;}
    [[nodiscard]] std::set<std::string> get_acceptors() const override { return acceptors;}
    [[nodiscard]] std::set<std::string> get_possible_links() const override { return possible_links;}
};

struct MANData: Data {
    std::set<std::string> donors = {"O2", "O3", "O6"};
    std::set<std::string> acceptors = {"C1"};
    std::set<std::string> possible_links = {"NAG", "MAN", };
    Sails::ResidueType type = Sails::ResidueType::NAG;
    [[nodiscard]] std::set<std::string> get_donors() const override { return donors;}
    [[nodiscard]] std::set<std::string> get_acceptors() const override { return acceptors;}
    [[nodiscard]] std::set<std::string> get_possible_links() const override { return possible_links;}
};

struct SailsData {
    SailsData() {
        auto* asndata = new ASNData();
        data_map.insert({Sails::ResidueType::ASN, asndata});

        auto* nagdata = new NAGData();
        data_map.insert({Sails::ResidueType::NAG, nagdata});

        auto* bmadata = new BMAData();
        data_map.insert({Sails::ResidueType::BMA, bmadata});

        auto* mandata = new MANData();
        data_map.insert({Sails::ResidueType::MAN, mandata});

    };

    std::map<Sails::ResidueType, Data*> data_map;

    ~SailsData() {
        for(auto& [fst, snd]: data_map) {
            delete snd;
        }
    };
};


#endif //SAILS_DATA_H
