/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH 1

#include "base_homomorphism_model.hh"

#include <vector>
#include <set>

class HomomorphismModel : public BaseHomomorphismModel
{
    protected:
        // Pairwise edge compatibilities between patern and target.
        std::vector<std::vector<bool> > edge_label_compatibility;

    public:
        HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params);
        auto target_edge_label(int t, int u) const -> int;
        using BaseHomomorphismModel::check_edge_label_compatibility;
        auto check_edge_label_compatibility(int p_lid, int t_lid) const -> bool;
};

#endif
