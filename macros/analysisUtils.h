#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

// Tiny helpers shared by step1.cc, step2.cc and buildtree.cc.
//
// Command-line convention:
//   ./program [inputFile] [--key=value ...]
// The input file is the first argument that does not start with "--".

#include <cstdlib>
#include <string>

// Value of --key=value, or `def` if the option was not given.
inline std::string getOption(int argc, char* argv[],
                             const std::string& key,
                             const std::string& def) {
    const std::string prefix = "--" + key + "=";
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s.rfind(prefix, 0) == 0)
            return s.substr(prefix.size());
    }
    return def;
}

// Numeric version of getOption.
inline double getNumber(int argc, char* argv[],
                        const std::string& key,
                        double def) {
    const std::string v = getOption(argc, argv, key, "");
    return v.empty() ? def : std::atof(v.c_str());
}

// First argument that is not an option.
inline std::string getInputFile(int argc, char* argv[],
                                const std::string& def) {
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s.rfind("--", 0) != 0)
            return s;
    }
    return def;
}

// data/d0_inclusive.root -> "_inclusive"
inline std::string tagFromFile(const std::string& inputName) {
    std::string stem = inputName.substr(inputName.find_last_of('/') + 1);
    if (stem.size() > 5 && stem.substr(stem.size() - 5) == ".root")
        stem.resize(stem.size() - 5);
    if (stem.rfind("d0_", 0) == 0)
        stem = stem.substr(3);
    return "_" + stem;
}

// Output-name tag: --tag=name if given, otherwise derived from the file name.
inline std::string getTag(int argc, char* argv[],
                          const std::string& inputName) {
    const std::string t = getOption(argc, argv, "tag", "");
    return t.empty() ? tagFromFile(inputName) : "_" + t;
}

#endif
