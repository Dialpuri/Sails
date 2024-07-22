//
// Created by Jordan Dialpuri on 20/07/2024.
//

#ifndef SAILS_CIF_H
#define SAILS_CIF_H

#include <string>
#include <utility>

namespace Sails {
    struct LinkRecord {
        LinkRecord(
            std::string id,
            gemmi::Chain &chain1,
            gemmi::Chain &chain2,
            gemmi::Residue &residue1,
            gemmi::Residue &residue2,
            gemmi::Atom &atom1,
            gemmi::Atom &atom2
        ): id(std::move(id)), chain1(chain1), chain2(chain2), residue1(residue1), residue2(residue2), atom1(atom1),
           atom2(atom2) {
        }

        static std::vector<std::string> tags() {
            return {
                "_struct_conn.id",
                "_struct_conn.conn_type_id",
                "_struct_conn.pdbx_leaving_atom_flag",
                "_struct_conn.pdbx_PDB_id",
                "_struct_conn.ptnr1_label_asym_id",
                "_struct_conn.ptnr1_label_comp_id",
                "_struct_conn.ptnr1_label_seq_id",
                "_struct_conn.ptnr1_label_atom_id",
                "_struct_conn.pdbx_ptnr1_label_alt_id",
                "_struct_conn.pdbx_ptnr1_PDB_ins_code",
                "_struct_conn.pdbx_ptnr1_standard_comp_id",
                "_struct_conn.ptnr1_symmetry",
                "_struct_conn.ptnr2_label_asym_id",
                "_struct_conn.ptnr2_label_comp_id",
                "_struct_conn.ptnr2_label_seq_id",
                "_struct_conn.ptnr2_label_atom_id",
                "_struct_conn.pdbx_ptnr2_label_alt_id",
                "_struct_conn.pdbx_ptnr2_PDB_ins_code",
                "_struct_conn.ptnr1_auth_asym_id",
                "_struct_conn.ptnr1_auth_comp_id",
                "_struct_conn.ptnr1_auth_seq_id",
                "_struct_conn.ptnr2_auth_asym_id",
                "_struct_conn.ptnr2_auth_comp_id",
                "_struct_conn.ptnr2_auth_seq_id",
                "_struct_conn.ptnr2_symmetry",
                "_struct_conn.pdbx_ptnr3_label_atom_id",
                "_struct_conn.pdbx_ptnr3_label_seq_id ",
                "_struct_conn.pdbx_ptnr3_label_comp_id",
                "_struct_conn.pdbx_ptnr3_label_asym_id",
                "_struct_conn.pdbx_ptnr3_label_alt_id",
                "_struct_conn.pdbx_ptnr3_PDB_ins_code",
                "_struct_conn.details",
                "_struct_conn.pdbx_dist_value",
                "_struct_conn.pdbx_value_order",
                "_struct_conn.pdbx_role",
            };
        }

        std::vector<std::string> labels() {
            double distance = (atom1.pos - atom2.pos).length();
            std::string res1_seqid = residue1.seqid.str();
            std::string res2_seqid = residue2.seqid.str();

            return {
                id,
                "covale",
                "one",
                "?",
                chain1.name,
                residue1.name,
                res1_seqid,
                atom1.name,
                "?",
                "?",
                "?",
                "1_555",
                chain2.name,
                residue2.name,
                ".",
                atom2.name,
                "?",
                "?",
                chain1.name,
                residue1.name,
                res1_seqid,
                chain2.name,
                residue2.name,
                res2_seqid,
                "1_555",
                "?",
                "?",
                "?",
                "?",
                "?",
                "?",
                "?",
                std::to_string(distance),
                "?",
                "N-Glycosylation"
            };
        }

        gemmi::Chain chain1;
        gemmi::Chain chain2;
        gemmi::Residue residue1;
        gemmi::Residue residue2;
        gemmi::Atom atom1;
        gemmi::Atom atom2;

        std::string id;
        std::string conn_type_id;
        std::string details;
        std::string pdbx_dist_value;
        std::string pdbx_role;
    };
}

#endif //SAILS_CIF_H
