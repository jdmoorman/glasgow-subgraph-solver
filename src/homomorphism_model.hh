/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH 1

#include "base_homomorphism_model.hh"
// #include "formats/input_graph.hh"
// #include "svo_bitset.hh"
// #include "homomorphism.hh"
// #include "homomorphism_domain.hh"
// #include "proof.hh"

// #include <memory>
#include <set>

class HomomorphismModel : public BaseHomomorphismModel
{
    public:
        using PatternAdjacencyBitsType = uint8_t;

        // const unsigned max_graphs;
        // unsigned pattern_size, target_size;
        //
        // auto has_less_thans() const -> bool;
        // std::vector<std::pair<unsigned, unsigned> > pattern_less_thans_in_convenient_order;

        HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params);
        ~HomomorphismModel();

        // auto pattern_vertex_for_proof(int v) const -> NamedVertex;
        // auto target_vertex_for_proof(int v) const -> NamedVertex;
        //sbool;
        //
        // auto pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType;
        // auto pattern_graph_row(int g, int p) const -> const SVOBitset &;
        // auto target_graph_row(int g, int t) const -> const SVOBitset &;
        //
        // auto forward_target_graph_row(int t) const -> const SVOBitset &;
        // auto reverse_target_graph_row(int t) const -> const SVOBitset &;
        //
        // auto pattern_degree(int g, int p) const -> unsigned;
        // auto target_degree(int g, int t) const -> unsigned;
        // auto largest_target_degree() const -> unsigned;

        auto has_vertex_labels() const -> bool;
        // TODO: Check, fix.
        auto has_edge_labels() const -> bool;
        // auto directed() const -> bool;
        // auto pattern_vertex_label(int p) const -> int;
        // auto target_vertex_label(int p) const -> int;
        auto pattern_edge_label(int p, int q) const -> int;
        auto target_edge_label(int t, int u) const -> int;
        auto edge_label_compatibility(int p_lid, int t_lid) const -> bool;
        // auto check_edge_label_compatibility(const std::multiset<std::string>&, const std::multiset<std::string>&) const -> bool;


        // auto pattern_has_loop(int p) const -> bool;
        // auto target_has_loop(int t) const -> bool;
        //
        // auto initialise_domains(std::vector<HomomorphismDomain> & domains) const -> bool;
};

#endif
