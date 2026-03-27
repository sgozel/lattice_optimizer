// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SG_UTILS_H
#define SG_UTILS_H

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <stdexcept>
#include "nlohmann/json.hpp"
using json = nlohmann::json;


// necessary for histogram optimization
static constexpr long long BASE = 85LL;  // must be strictly > number of bonds


struct Site {
    int    index;
    double x, y, z;
};
std::vector<Site> sites_list;


struct Bond {
	int first;
	int second;
	std::string coupling;
};

// -----------------------------------------------------------------------------
// Basic metrics
// perm[v] = position of vertex v
// -----------------------------------------------------------------------------

int bandwidth(const std::vector<int>& perm, const std::vector<Bond>& bonds)
{
    int bw = 0;
    for (size_t k = 0; k < bonds.size(); ++k) { 
		const int i = bonds[k].first;
		const int j = bonds[k].second;
        bw = std::max(bw, std::abs(perm[i] - perm[j]));
    }
    return bw;
}


int total_dist(const std::vector<int>& perm, const std::vector<Bond>& bonds)
{
	int dist = 0;
    for (int k = 0; k < (int)bonds.size(); ++k) {
		const int i = bonds[k].first;
		const int j = bonds[k].second;
        dist += std::abs(perm[i] - perm[j]);
    }
    return dist;
}


int total_ops(const std::vector<int>& perm, const std::vector<Bond>& bonds)
{
    int ops = 0;
    for (size_t k = 0; k < bonds.size(); ++k) {
		const int i = bonds[k].first;
		const int j = bonds[k].second;
        ops += 2 * std::abs(perm[i] - perm[j]) - 1;
    }
    return ops;
}

// -----------------------------------------------------------------------------
// Histogram
// -----------------------------------------------------------------------------

using Histogram = std::vector<int>;

Histogram compute_histogram(const std::vector<int>& perm, const std::vector<Bond>& bonds)
{
	unsigned int n = perm.size();
    Histogram hist(n, 0);
    for (size_t k = 0; k < bonds.size(); ++k) {
		const int i = bonds[k].first;
		const int j = bonds[k].second;
        hist[std::abs(perm[i] - perm[j])]++;
    }
    return hist;
}


// -----------------------------------------------------------------------------
// Greedy BFS warm-start ordering
// -----------------------------------------------------------------------------

std::vector<int> greedy_ordering(const std::vector<Bond>& bonds, const unsigned int n)
{
    std::vector<std::vector<int>> adj(n);
    for (size_t k = 0; k < bonds.size(); ++k) {
		const int i = bonds[k].first;
		const int j = bonds[k].second;
        adj[i].push_back(j);
        adj[j].push_back(i);
    }
    int start = (int)(std::max_element(adj.begin(), adj.end(),
        [](const auto& a, const auto& b){ return a.size() < b.size(); })
        - adj.begin());

    std::vector<int> perm(n, 0);
    std::vector<bool> visited(n, false);
    std::vector<int> queue = {start};
    visited[start] = true;
    int pos = 0;
    while (!queue.empty()) {
        std::vector<int> next;
        for (int u : queue) {
            perm[u] = pos++;
            for (int w : adj[u]) {
                if (!visited[w]) {
					visited[w] = true;
					next.push_back(w);
				}
			}
        }
        queue = next;
    }
    return perm;
}

// -----------------------------------------------------------------------------
// Simulated Annealing  (objective returns long long)
// -----------------------------------------------------------------------------

struct SAResult {
    std::vector<int> perm;
    long long score;
};

struct SAParams {
	int n_iter;
    double T_start;
    double T_end;
    int seed;
};


template<typename ObjFn>
SAResult simulated_annealing(
    ObjFn objective_fn,
    const std::vector<int>& init_perm,
    const SAParams& saparams)
{
	const unsigned int n = init_perm.size();
	
    std::mt19937 rng(saparams.seed);
    std::uniform_int_distribution<int> pos_dist(0, n-1);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);

    // ordering[pos] = vertex — swaps are O(1) and trivially undoable
    std::vector<int> ordering(n);
    for (unsigned int v = 0; v < n; ++v) {
		ordering[init_perm[v]] = v;
	}

    std::vector<int> best_ordering = ordering;
    long long current_score = objective_fn(init_perm);
    long long best_score = current_score;

    const double log_ratio = std::log(saparams.T_end / saparams.T_start);

    for (int it = 0; it < saparams.n_iter; ++it)
    {
        const double T = saparams.T_start * std::exp(log_ratio * it / saparams.n_iter);

        int a = pos_dist(rng);
        int b = pos_dist(rng);
        if (a == b) {
			continue;
		}

        std::swap(ordering[a], ordering[b]);

        std::vector<int> new_perm(n);
        for (int p = 0; p < n; ++p) {
			new_perm[ordering[p]] = p;
		}

        const long long new_score = objective_fn(new_perm);
        const double delta = (double)(new_score - current_score);

        if (delta < 0.0 || real_dist(rng) < std::exp(-delta / T)) {
            current_score = new_score;
            if (current_score < best_score) {
                best_score = current_score;
                best_ordering = ordering;
            }
        } else {
            std::swap(ordering[a], ordering[b]);  // undo
        }
    }

    std::vector<int> best_perm(n);
    for (int p = 0; p < n; ++p) {
		best_perm[best_ordering[p]] = p;
	}
    return {best_perm, best_score};
}


// -----------------------------------------------------------------------------
// Output helpers
// -----------------------------------------------------------------------------


void print_perm(const std::string& label, 
				const std::vector<Bond>& bonds, 
				const std::vector<int>& perm)
{
	const unsigned int n = perm.size();
	
    if (!label.empty()) {
		std::cout << label << std::endl;
	}
    std::cout << "  vertex->pos : [";
    for (unsigned int v = 0; v < n; ++v) {
		std::cout << perm[v] << (v<n-1?",":"]");
	}
	std::cout << std::endl;
    std::vector<int> ord(n);
    for (unsigned int v = 0; v < n; ++v) {
		ord[perm[v]] = v;
	}
    std::cout << "  pos->vertex : [";
    for (unsigned int p = 0; p < n; ++p) {
		std::cout << ord[p] << (p<n-1?",":"]");
	}
	std::cout << std::endl;
    
    const int bw = bandwidth(perm, bonds);
    const int totops = total_ops(perm, bonds);
    const Histogram hist = compute_histogram(perm, bonds);
    
    std::cout << "  Bandwidth : " << bw << std::endl;
    std::cout << "  Total ops : " << totops << std::endl;
    std::cout << "  Histogram (d : count : ops_size)" << std::endl;
    for (unsigned int d = bw; d >= 1; --d) {
        if (hist[d] > 0) {
            std::cout << "    d=" << std::setw(2) << d
                      << ", count=" << std::setw(3) << hist[d]
                      << ", ops_size=" << (2*d-1) << std::endl;
		}
	}
}


void print_bw_bond_starts(const std::vector<Bond>& bonds, 
						  const std::vector<int>& perm, 
						  int target_bw)
{
    std::vector<int> starts;
    for (size_t k = 0; k < bonds.size(); ++k) {
        const int i = bonds[k].first;
        const int j = bonds[k].second;
        if (std::abs(perm[i] - perm[j]) == target_bw) {
            starts.push_back(std::min(perm[i], perm[j]));
		}
    }
    std::sort(starts.begin(), starts.end());
    std::cout << "  BW-bond starts (sorted): [";
    for (size_t k = 0; k < starts.size(); ++k) {
        std::cout << starts[k] << (k < (starts.size()-1) ? "," : "]");
	}
    std::cout << std::endl;
}


void print_table(const std::vector<Bond>& bonds, const std::vector<int>& perm)
{
    std::cout << "\n";
    std::cout << std::setw(4)  << "k"
              << std::setw(10) << "original"
              << std::setw(8)  << "|i-j|"
              << std::setw(6)  << "ops"
              << std::setw(12) << "transformed"
              << std::setw(10) << "|i-j|'"
              << std::setw(8)  << "ops'"
              << "\n";
    std::cout << std::string(58, '-') << "\n";
    for (size_t k = 0; k < bonds.size(); ++k) {
        const int i = bonds[k].first;
        const int j = bonds[k].second;
        const int orig_d = std::abs(i - j);
        int ni = perm[i];
        int nj = perm[j];
        if (ni > nj) {
			std::swap(ni, nj);
		}
        const int new_d = nj - ni;
        
        std::string os = "("+std::to_string(i)+","+std::to_string(j)+")";
        std::string ns = "("+std::to_string(ni)+","+std::to_string(nj)+")";
        std::cout << std::setw(4)  << k
                  << std::setw(10) << os
                  << std::setw(8)  << orig_d
                  << std::setw(6)  << (2*orig_d-1)
                  << std::setw(12) << ns
                  << std::setw(10) << new_d
                  << std::setw(8)  << (2*new_d-1)
                  << "\n";
    }
}

// -----------------------------------------------------------------------------
// Print transformed bond file (same format as input bond file)
// -----------------------------------------------------------------------------

void print_transformed_bond_file(const std::vector<Bond>& bonds, 
								 const std::vector<Site>& sites_list, 
								 const std::vector<int>& perm)
{
	const unsigned int n = perm.size();
	
    // Build inverse permutation: pos -> original vertex index
    std::vector<int> inv_perm(n);
    for (unsigned int v = 0; v < n; ++v) {
		inv_perm[perm[v]] = (int)v;
	}

    std::cout << "[sites]=" << n << std::endl;
    for (unsigned int p = 0; p < n; ++p) {
        const Site& s = sites_list[inv_perm[p]];
        std::cout << p << " " << s.x << " " << s.y << " " << s.z << std::endl;
    }
    std::cout << "[interactions]=" << bonds.size() << std::endl;
    for (size_t k = 0; k < bonds.size(); ++k) {
        const int i = bonds[k].first;
        const int j = bonds[k].second;
        const int ni = perm[i];
        const int nj = perm[j];
        std::cout << bonds[k].coupling 
				  << " (" << std::min(ni, nj) << ", " << std::max(ni, nj) 
				  << ")" << std::endl;
    }
}

// -----------------------------------------------------------------------------
// Input parsing
// -----------------------------------------------------------------------------

void load_bonds_from_file(const std::string& filename, 
						  std::vector<Site>& sites_list, 
						  std::vector<Bond>& bonds)
{
    std::ifstream f(filename);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open bonds file: " + filename);
	}

    std::string line;
    unsigned int n;

    // First line: [sites]=<value>
    if (!std::getline(f, line)) {
        throw std::runtime_error("Bonds file is empty: " + filename);
	}
    {
        auto eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error("Invalid format for sites line: " + line);
		}
        n = std::stoi(line.substr(eq + 1));
    }

    // Next n lines: site definitions
    // Format: <index> <x> <y> <z>
    sites_list.clear();
    sites_list.reserve(n);
    for (int s = 0; s < n; ++s) {
        if (!std::getline(f, line)) {
            throw std::runtime_error("Unexpected end of file while reading sites.");
		}
        std::istringstream iss(line);
        Site site;
        if (!(iss >> site.index >> site.x >> site.y >> site.z)) {
            throw std::runtime_error("Invalid site line: " + line);
		}
        sites_list.push_back(site);
    }
    if ((int)sites_list.size() != n) {
        throw std::runtime_error("Expected " + std::to_string(n) +
                                 " sites but read " + std::to_string(sites_list.size()));
	}

    // Next line: [interactions]=<value>
    if (!std::getline(f, line)) {
        throw std::runtime_error("Bonds file missing interactions count line.");
	}
    int num_bonds = 0;
    {
        auto eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error("Invalid format for interactions line: " + line);
		}
        num_bonds = std::stoi(line.substr(eq + 1));
    }

    // Remaining lines: couplingName (i, j)
    bonds.clear();
	bonds.reserve(num_bonds);
	while (std::getline(f, line)) {
		if (line.empty()) {
			continue;
		}
		auto open  = line.find('(');
		auto comma = line.find(',');
		auto close = line.find(')');
		if (open == std::string::npos || comma == std::string::npos || close == std::string::npos) {
			throw std::runtime_error("Invalid bond line: " + line);
		}

		// Extract coupling name: everything before '(' stripped of whitespace
		std::string coupling = line.substr(0, open);
		coupling.erase(std::remove_if(coupling.begin(), coupling.end(), ::isspace), coupling.end());
		if (coupling.empty()) {
			throw std::runtime_error("Missing coupling name in bond line: " + line);
		}

		int i = std::stoi(line.substr(open + 1, comma - open - 1));
		int j = std::stoi(line.substr(comma + 1, close - comma - 1));
		bonds.emplace_back(Bond{i, j, coupling});
	}
    
    if ((int)bonds.size() != num_bonds) {
        throw std::runtime_error("Expected " + std::to_string(num_bonds) +
                                 " bonds but read " + std::to_string(bonds.size()));
	}
	/*
	if (BASE<=bonds.size()) {
		throw std::runtime_error("BASE<=bonds.size(). Increase BASE.");
	}*/
}

// -----------------------------------------------------------------------------
// Read input .json file
// -----------------------------------------------------------------------------

json read_input_file(int argc, char* argv[])
{
	std::ifstream inputFileStream;
	json inputParam;
    try {
        if (argc == 2) {
            inputFileStream.open(argv[1]);
            if (!inputFileStream.is_open()) {
                throw std::runtime_error("Cannot open input file: " + std::string(argv[1]));
            }
            inputFileStream >> inputParam;
            inputFileStream.close();
        } else {
            throw std::runtime_error("No input file provided.");
        }
    }
    catch (...) {
        std::cerr << "Caught exception at top level." << std::endl;
        std::abort();
    }
    return inputParam;
}

void read_json_file(const json& inputParam, 
					std::vector<Site>& sites_list,
					std::vector<Bond>& bonds)
{
	// Read bonds filename from JSON and load bonds
    if (!inputParam.contains("bonds_filename")) {
        throw std::runtime_error("Missing bonds_filename in input .json file.");
	}
    std::string bonds_filename = inputParam["bonds_filename"];
    
    load_bonds_from_file(bonds_filename, sites_list, bonds);
}

#endif
