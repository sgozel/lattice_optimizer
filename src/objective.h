// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SG_OBJECTIVE_H
#define SG_OBJECTIVE_H

#include <vector>

#include "utils.h"

// -----------------------------------------------------------------------------
// Total ops objective
// -----------------------------------------------------------------------------

long long constrained_total_ops_score(const std::vector<int>& perm, 
									  const std::vector<Bond>& bonds, 
									  int target_bw)
{
    int bw = bandwidth(perm, bonds);
    if (bw > target_bw) {
        return 1000000000000LL + (long long)bw * 1000000000LL;
	}
    return (long long)total_ops(perm, bonds);
}


// -----------------------------------------------------------------------------
// Bandwidth objective
// -----------------------------------------------------------------------------

long long constrained_bandwidth_score(const std::vector<int>& perm, 
									  const std::vector<Bond>& bonds,
									  int target_totops)
{
	int to = total_ops(perm, bonds);
    if (to > target_totops) {
        return 1000000000000LL + (long long)to * 1000000000LL;
	}
    return (long long)bandwidth(perm, bonds);
}


// -----------------------------------------------------------------------------
// Lexicographic histogram objective
//
// hist[d] = number of bonds with |pi[i]-pi[j]| == d
//
// We want to minimise (hist[BW], hist[BW-1], ..., hist[1]) lexicographically.
//
// Scalar encoding for SA:
//   score = sum_{d=1}^{BW} hist[d] * BASE^d
// where BASE > total number of bonds.
// This ensures that reducing hist[d] by 1 is always better than zeroing
// all hist[d'] for every d' < d.
// With 40 bonds, BASE = 41 suffices.
//
// If BW > target_bw, return a large penalty.
// -----------------------------------------------------------------------------

long long histogram_score(const Histogram& hist, int target_bw)
{
    long long score  = 0;
    long long weight = BASE;   // BASE^1 for d=1
    for (int d = 1; d <= target_bw; ++d) {
        score  += (long long)hist[d] * weight;
        weight *= BASE;
    }
    return score;
}

long long constrained_histogram_score(const std::vector<int>& perm, 
									  const std::vector<Bond>& bonds,
									  int target_bw)
{
    int bw = bandwidth(perm, bonds);
    if (bw > target_bw) {
        return 1000000000000LL + (long long)bw * 1000000000LL;
	}
    return histogram_score(compute_histogram(perm, bonds), target_bw);
}


// -----------------------------------------------------------------------------
// Objective: minimize starting indices of BW-distance bonds
//
// Subject to:
//   - bandwidth(perm) <= target_bw
//   - hist[target_bw]   == target_hist_bw    (count at BW preserved)
//   - hist[target_bw-1] == target_hist_bw_m1 (count at BW-1 preserved)
//
// For all bonds at distance target_bw, collect their starting index
// (= min of the two endpoint positions), sort descending, then encode
// as a weighted sum so that minimizing the score pushes the largest
// starting index down first, then the second largest, etc.
//
// Encoding: score = sum_k starts_sorted_desc[k] * N^k
// where N > n (so each "digit" can hold any position index).
// We use N = n+1.
// -----------------------------------------------------------------------------

long long starting_index_score(
    const std::vector<int>& perm,
    const std::vector<Bond>& bonds,
    int target_bw,
    int target_hist_bw,
    int target_hist_bw_m1)
{
    // Bandwidth penalty
    int bw = bandwidth(perm, bonds);
    if (bw > target_bw) {
        return 2000000000000LL + (long long)bw * 1000000000LL;
	}

    // Histogram constraints for BW and BW-1
    Histogram hist = compute_histogram(perm, bonds);
    if (hist[target_bw] != target_hist_bw) {
        return 1500000000000LL + (long long)std::abs(hist[target_bw] - target_hist_bw) * 1000000000LL;
	}
    if (target_bw >= 2 && hist[target_bw - 1] != target_hist_bw_m1) {
        return 1000000000000LL + (long long)std::abs(hist[target_bw - 1] - target_hist_bw_m1) * 1000000000LL;
	}

    // Collect starting indices of BW-distance bonds
    std::vector<int> starts;
    starts.reserve(target_hist_bw);
    for (size_t k = 0; k < bonds.size(); ++k) {
        const int i = bonds[k].first;
        const int j = bonds[k].second;
        const int d = std::abs(perm[i] - perm[j]);
        if (d == target_bw) {
            starts.push_back(std::min(perm[i], perm[j]));
		}
    }

    // Sort descending: we minimize the largest starting index first
    std::sort(starts.begin(), starts.end(), std::greater<int>());

    // Encode as positional score: starts[0]*N^0 + starts[1]*N^1 + ...
    // (starts[0] is the largest, so it has the highest weight in the lex order)
    long long score = 0;
    long long weight = 1;
    const long long N = (long long)(perm.size() + 1);
    for (int s : starts) {
        score += (long long)s * weight;
        weight *= N;
    }
    return score;
}

#endif
