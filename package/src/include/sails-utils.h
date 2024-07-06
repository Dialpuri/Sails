//
// Created by jordan on 05/07/24.
//

#ifndef SAILS_SAILS_UTILS_H
#define SAILS_SAILS_UTILS_H

#include <iostream>
#include <filesystem>

#include "gemmi/math.hpp"

namespace Sails::Utils {
    inline bool check_file_exists(const std::string& path) {
        return std::filesystem::exists(path);
    }

    inline void print_vector(const gemmi::Vec3& vec) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }

} // namespace Sails::Utils


#endif //SAILS_SAILS_UTILS_H
