#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <algorithm>
#include <iostream>

namespace pdf {

class FastInterpolator {
public:
    FastInterpolator() = default;

    // Precompute interpolation weights mapping source_grid -> target_grid
    void Initialize(const std::vector<double>& source_grid, const std::vector<double>& target_grid) {
        indices_.resize(target_grid.size());
        weights_.resize(target_grid.size());

        for (size_t i = 0; i < target_grid.size(); ++i) {
            double x = target_grid[i];
            
            // Find position in source grid
            auto it = std::lower_bound(source_grid.begin(), source_grid.end(), x);
            
            if (it == source_grid.begin()) {
                indices_[i] = 0;
                weights_[i] = 0.0; // Use source[0]
            } else if (it == source_grid.end()) {
                indices_[i] = source_grid.size() - 2;
                weights_[i] = 1.0; // Use source[last]
            } else {
                size_t idx = std::distance(source_grid.begin(), it) - 1;
                indices_[i] = idx;
                
                double x1 = source_grid[idx];
                double x2 = source_grid[idx+1];
                weights_[i] = (x - x1) / (x2 - x1);
            }
        }
    }

    // Interpolate values from source_y to target_y using precomputed weights
    void Interpolate(const std::vector<double>& source_y, std::vector<double>& target_y) const {
        if (target_y.size() != indices_.size()) {
            target_y.resize(indices_.size());
        }
        
        for (size_t i = 0; i < indices_.size(); ++i) {
            size_t idx = indices_[i];
            double w = weights_[i];
            target_y[i] = source_y[idx] * (1.0 - w) + source_y[idx + 1] * w;
        }
    }

private:
    std::vector<size_t> indices_;
    std::vector<double> weights_;
};

} // namespace pdf

#endif // INTERPOLATION_HPP
