/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "IO.h++"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace LED::IO {
bool Output::create_directory_if_needed(const std::string& filepath) const {
    try {
        std::filesystem::path path(filepath);
        auto parent_path = path.parent_path();
        if (!parent_path.empty() && !std::filesystem::exists(parent_path)) {
            return std::filesystem::create_directories(parent_path);
        }
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
        return false;
    }
}

void Output::InitOutput(const std::string& filename, const std::string& headline) {
    fOutputFilePath = filename;

    if (fOutputFilePath.empty()) {
        std::cout << headline;
    } else {
        if (!create_directory_if_needed(fOutputFilePath)) {
            std::cerr << "Failed to create directory for: " << fOutputFilePath << std::endl;
            fOutputFilePath.clear();
            return;
        }

        std::ofstream f(fOutputFilePath);
        if (!f.is_open()) {
            std::cerr << "File cannot be opened: " << fOutputFilePath << std::endl;
            fOutputFilePath.clear();
        } else {
            f << headline;
        }
    }
}

void Output::AddToOutput(const double n1, const double n2, const double n3) {
    if (fOutputFilePath.empty()) {
        std::cout << n1 << " " << n2 << " " << n3 << std::endl;
    } else {
        std::ofstream f(fOutputFilePath, std::ios::app);
        if (!f.is_open()) {
            std::cerr << "File cannot be opened: " << fOutputFilePath << std::endl;
            fOutputFilePath.clear();
        } else {
            f << n1 << " " << n2 << " " << n3 << std::endl;
        }
    }
}

void Output::AddToOutput2(const double n1, const double n2) {
    if (fOutputFilePath.empty()) {
        std::cout << n1 << " " << n2 << std::endl;
    } else {
        std::ofstream f(fOutputFilePath, std::ios::app);
        if (!f.is_open()) {
            std::cerr << "File cannot be opened: " << fOutputFilePath << std::endl;
            fOutputFilePath.clear();
        } else {
            f << n1 << " " << n2 << std::endl;
        }
    }
}
} // namespace LED::IO
