// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <vector>
#include <omp.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "utils.h"
#include "objective.h"

int main(int argc, char* argv[])
{
	json inputParam = read_input_file(argc, argv);
	std::vector<Site> sites_list;
	std::vector<Bond> bonds;
	read_json_file(inputParam, sites_list, bonds);
	const unsigned int n = sites_list.size();
	
    std::cout << "Loaded n=" << n 
			  << ", bonds=" << bonds.size()
              << " from '" << inputParam["bonds_filename"] 
              << std::endl;
	
    auto t_start = std::chrono::high_resolution_clock::now();
	
	// ---- Tunable parameters ----
    int pass1_seeds = 400;
    SAParams saparams_pass1_global;
    saparams_pass1_global.n_iter  = 500000;
    saparams_pass1_global.T_start = 8.0;
    saparams_pass1_global.T_end   = 0.001;
    
    int pass2_seeds  = 400;
    SAParams saparams_pass2_global;
    saparams_pass2_global.n_iter  = 500000;
    saparams_pass2_global.T_start = 20.0;
    saparams_pass2_global.T_end   = 0.002;
    // ----------------------------
	
	// Baseline - identity permutation - starting indexing
    std::vector<int> identity(n);
    std::iota(identity.begin(), identity.end(), 0);
    print_perm("=== ORIGINAL (identity) ===", bonds, identity);
	
	//------------------------------------------------------------------
	// --- Pass 1: minimize total ops ---
	//------------------------------------------------------------------
    std::cout << "\n=== PASS 1: minimize total ops ===" << std::endl;

	std::cout << "max_threads: " << omp_get_max_threads() << std::endl;
	
    int best_totops = std::numeric_limits<int>::max();
    std::vector<std::vector<int>> totops_pool;
    const auto init = greedy_ordering(bonds, n);

	unsigned int cpt_seed = 0;

	#pragma omp parallel for schedule(static)
    for (int seed = 0; seed < pass1_seeds; ++seed) {
		if ((cpt_seed > 0) && (cpt_seed%10 == 0)) {
			std::cout << "cpt_seed = " << cpt_seed << std::endl;
		}
		SAParams saparams_pass1 = saparams_pass1_global;
        saparams_pass1.seed = seed;
        auto result = simulated_annealing(
            [&bonds](const std::vector<int>& p) -> long long { return total_ops(p, bonds); },
            init, 
            saparams_pass1);

        int totops = total_ops(result.perm, bonds);
        if (totops < best_totops) {
			#pragma omp critical
			{
				best_totops = totops;
				totops_pool.clear();
				std::cout << "  seed=" << std::setw(3) << seed
						  << "  new best totops=" << best_totops
						  << "  total_ops=" << total_ops(result.perm, bonds) 
						  << std::endl;
			}
        }
        if (totops == best_totops) {
			#pragma omp critical
			{
				totops_pool.push_back(result.perm);
			}
		}
		#pragma omp atomic
		cpt_seed += 1;
    }
    
    std::cout << "  Best tot ops=" << best_totops
              << "  Pool size=" << totops_pool.size() << std::endl;
	
	//------------------------------------------------------------------
	// --- Pass 2: bandwidth minimization ---
	//------------------------------------------------------------------
	
	const int target_totops = best_totops;
    std::cout << "\n=== PASS 2: bandwidth minimization"
              << " subject to tot ops<=" << target_totops << " ===" << std::endl;

    long long best_score2 = std::numeric_limits<long long>::max();
    std::vector<int> best_perm2 = totops_pool[0];

	cpt_seed = 0;
	
	#pragma omp parallel for schedule(static)
    for (int seed = 0; seed < pass2_seeds; ++seed) {
		if ((cpt_seed > 0) && (cpt_seed%10 == 0)) {
			std::cout << "cpt_seed = " << cpt_seed << std::endl;
		}
		SAParams saparams_pass2 = saparams_pass2_global;
		saparams_pass2.seed = seed + 2000;
        auto& init2 = totops_pool[seed % totops_pool.size()];
        auto result = simulated_annealing(
            [&bonds,target_totops](const std::vector<int>& p) -> long long {
				//======================================
				// HERE CHOICE OF OBJECTIVE FUNCTION
				//======================================
                return constrained_bandwidth_score(p, bonds, target_totops); // bandwidth
            },
            init2, 
            saparams_pass2);

        if (total_ops(result.perm, bonds) <= target_totops && result.score < best_score2) {
            #pragma omp critical
            {
				best_score2 = result.score;
				best_perm2  = result.perm;
				Histogram hist = compute_histogram(best_perm2, bonds);
				std::cout << "  seed=" << std::setw(3) << seed
						  << "  BW=" << bandwidth(best_perm2, bonds)
						  << "  total_ops=" << total_ops(best_perm2, bonds)
						  << "  hist=[";
				for (int d = bandwidth(best_perm2, bonds); d >= 1; --d) {
					std::cout << "d" << d << ":" << hist[d] << (d>1?" ":"]");
				}
				std::cout << std::endl;
			}
        }
        #pragma omp atomic
        cpt_seed += 1;
    }
	
	//------------------------------------------------------------------
	// Final result
	//------------------------------------------------------------------
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "\n=== FINAL RESULT (elapsed: "
              << std::fixed << std::setprecision(2) 
              << elapsed << "s) ===" << std::endl;
    print_perm("", bonds, best_perm2);
    print_table(bonds, best_perm2);
	
    std::cout << "\n=== TRANSFORMED BOND FILE ===" << std::endl;
    print_transformed_bond_file(bonds, sites_list, best_perm2);
	
	return 0;
}
