/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CROSSWORD_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CROSSWORD_HOMOMORPHISM_MODEL_HH 1

#include "base_homomorphism_model.hh"
#include "formats/input_graph.hh"

#include <map>
#include <string>
#include <set>

class CrosswordHomomorphismModel : public BaseHomomorphismModel
{
    protected:
        std::vector<std::string> target_vertex_names;
        std::vector<std::multiset<std::string>> id_to_pattern_edge_labels;

    public:
        CrosswordHomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params);
        using BaseHomomorphismModel::check_edge_label_compatibility;
        auto check_edge_label_compatibility(const int t_v1, const int t_v2, const int p_lid) const -> bool;
};

#endif
