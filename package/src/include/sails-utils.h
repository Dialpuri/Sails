//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_UTILS_H
#define SAILS_SAILS_UTILS_H

#include <iostream>
#include <filesystem>


namespace Sails::Utils {
    inline bool check_file_exists(const std::string& path) {
        return std::filesystem::exists(path);
    }
} // namespace Sails::Utils


#endif //SAILS_SAILS_UTILS_H
