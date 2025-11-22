#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string Trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return std::string();
    }
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::vector<std::string> SplitWhitespace(const std::string& text)
{
    std::istringstream in(text);
    std::vector<std::string> parts;
    std::string token;
    while (in >> token) {
        parts.push_back(token);
    }
    return parts;
}

struct AliasInfo {
    std::vector<std::string> volumes;
    std::map<std::string, std::string> metadata;
};

AliasInfo ParseAliasFile(const std::string& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Unable to open alias file: " + path);
    }

    AliasInfo info;
    std::string line;
    size_t line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        const auto comment = line.find('#');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }

        line = Trim(line);
        if (line.empty()) {
            continue;
        }

        const auto key_end = line.find_first_of(" \t");
        if (key_end == std::string::npos) {
            throw std::runtime_error("Line " + std::to_string(line_number) +
                                     " is missing a value: " + line);
        }

        const std::string key = line.substr(0, key_end);
        const std::string value = Trim(line.substr(key_end + 1));
        if (value.empty()) {
            throw std::runtime_error("Line " + std::to_string(line_number) +
                                     " has an empty value: " + line);
        }

        if (key == "DBLIST") {
            const auto volumes = SplitWhitespace(value);
            if (volumes.empty()) {
                throw std::runtime_error("DBLIST on line " + std::to_string(line_number) +
                                         " does not list any volumes");
            }
            info.volumes.insert(info.volumes.end(), volumes.begin(), volumes.end());
            continue;
        }

        const auto duplicate = info.metadata.find(key);
        if (duplicate != info.metadata.end()) {
            throw std::runtime_error("Duplicate key '" + key + "' on line " +
                                     std::to_string(line_number));
        }

        info.metadata[key] = value;
    }

    return info;
}

void PrintAliasInfo(const AliasInfo& info)
{
    std::cout << "Volumes (DBLIST):" << std::endl;
    if (info.volumes.empty()) {
        std::cout << "  <none>" << std::endl;
    } else {
        for (const auto& name : info.volumes) {
            std::cout << "  - " << name << std::endl;
        }
    }

    std::cout << std::endl << "Additional data:" << std::endl;
    if (info.metadata.empty()) {
        std::cout << "  <none>" << std::endl;
    } else {
        for (const auto& kv : info.metadata) {
            std::cout << "  " << kv.first << ": " << kv.second << std::endl;
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <alias-file.pal>" << std::endl;
        return 1;
    }

    try {
        const AliasInfo info = ParseAliasFile(argv[1]);
        PrintAliasInfo(info);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

