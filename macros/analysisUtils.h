#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

// Tiny helpers shared by sanity.cc, mass.cc, features.cc, buildtree.cc and
// the generators (jets-charm.cc, jets-inc.cc).
//
// Command-line convention:
//   ./program [inputFile] [--key=value ...]
// The input file is the first argument that does not start with "--".

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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

// index-th (0-based) argument that is not an option, or `def`.
inline std::string getPositional(int argc, char* argv[], int index,
                                 const std::string& def) {
    int count = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s.rfind("--", 0) == 0)
            continue;
        if (count == index)
            return s;
        ++count;
    }
    return def;
}

// data/d0_inclusive.root         -> "_inclusive"
// data/candidates_charm.root     -> "_charm"
inline std::string tagFromFile(const std::string& inputName) {
    std::string stem = inputName.substr(inputName.find_last_of('/') + 1);
    if (stem.size() > 5 && stem.substr(stem.size() - 5) == ".root")
        stem.resize(stem.size() - 5);
    if (stem.rfind("d0_", 0) == 0)
        stem = stem.substr(3);
    if (stem.rfind("candidates_", 0) == 0)
        stem = stem.substr(11);
    return "_" + stem;
}

// Output-name tag: --tag=name if given, otherwise derived from the file name.
inline std::string getTag(int argc, char* argv[],
                          const std::string& inputName) {
    const std::string t = getOption(argc, argv, "tag", "");
    return t.empty() ? tagFromFile(inputName) : "_" + t;
}

// Reject unknown or malformed options, so a typo such as --smear-pT=0.03 is
// an error instead of being silently ignored. Returns false (after printing
// what is wrong) if any "--" argument is not written --key=value or its key
// is not in `known`.
inline bool checkOptions(int argc, char* argv[],
                         const std::vector<std::string>& known) {
    bool ok = true;
    for (int i = 1; i < argc; ++i) {
        const std::string s = argv[i];
        if (s.rfind("--", 0) != 0)
            continue;
        const size_t eq = s.find('=');
        if (eq == std::string::npos) {
            std::cerr << "Option '" << s << "' must be written --key=value\n";
            ok = false;
            continue;
        }
        const std::string key = s.substr(2, eq - 2);
        if (std::find(known.begin(), known.end(), key) == known.end()) {
            std::cerr << "Unknown option --" << key << "\nKnown options:";
            for (const std::string& k : known)
                std::cerr << " --" << k;
            std::cerr << "\n";
            ok = false;
        }
    }
    return ok;
}

// One-line-per-run CSV. Replaces the row whose first `nKey` comma-separated
// fields match `row`, or appends it. If the file exists with a different
// header (older column layout) it is renamed to <path>.old first.
inline void upsertCsvRow(const std::string& path, const std::string& header,
                         const std::string& row, int nKey = 2) {
    const auto keyOf = [nKey](const std::string& l) {
        size_t pos = 0;
        for (int k = 0; k < nKey; ++k) {
            pos = l.find(',', pos);
            if (pos == std::string::npos)
                return l;
            ++pos;
        }
        return l.substr(0, pos);
    };

    std::vector<std::string> lines;
    bool headerOk = true;
    {
        std::ifstream in(path);
        std::string line;
        bool first = true;
        while (std::getline(in, line)) {
            if (first) {
                first = false;
                if (line != header) {
                    headerOk = false;
                    break;
                }
                continue;
            }
            if (!line.empty())
                lines.push_back(line);
        }
        if (!first && !headerOk) {
            in.close();
            std::rename(path.c_str(), (path + ".old").c_str());
            lines.clear();
        }
    }

    const std::string key = keyOf(row);
    std::ofstream out(path);
    out << header << "\n";
    bool replaced = false;
    for (const std::string& l : lines) {
        if (!replaced && keyOf(l) == key) {
            out << row << "\n";
            replaced = true;
        } else if (keyOf(l) != key) {
            out << l << "\n";
        }
    }
    if (!replaced)
        out << row << "\n";
}

#endif
