/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH 1

#include "formats/input_graph.hh"
#include "svo_bitset.hh"
#include "homomorphism.hh"
#include "homomorphism_domain.hh"
#include "proof.hh"

#include <memory>
#include <set>

class HomomorphismModel
{
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

        auto _build_exact_path_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
                unsigned number_of_exact_path_graphs, bool directed) -> void;

        auto _build_distance3_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void;

        auto _build_k4_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void;

        auto _check_degree_compatibility(
                int p,
                int t,
                unsigned graphs_to_consider,
                std::vector<std::vector<std::vector<int> > > & patterns_ndss,
                std::vector<std::vector<std::optional<std::vector<int> > > > & targets_ndss,
                bool do_not_do_nds_yet
                ) const -> bool;

        auto _check_edge_label_compatibility(const std::multiset<std::string>&, const std::multiset<std::string>&) const -> bool;

        auto _check_loop_compatibility(int p, int t) const -> bool;

        auto _check_label_compatibility(int p, int t) const -> bool;

        auto _label_counts_from_multiset(const std::multiset<std::string>&) const -> std::map<std::string, int>;

        auto _record_edge_labels(std::map<std::multiset<std::string>, int>& label_map, const InputGraph & graph, std::vector<int>& graph_edge_labels) -> void;

    public:
        using PatternAdjacencyBitsType = uint8_t;

        const unsigned max_graphs;
        unsigned pattern_size, target_size;

        auto has_less_thans() const -> bool;
        std::vector<std::pair<unsigned, unsigned> > pattern_less_thans_in_convenient_order;

        HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params);
        ~HomomorphismModel();

        auto pattern_vertex_for_proof(int v) const -> NamedVertex;
        auto target_vertex_for_proof(int v) const -> NamedVertex;

        auto prepare() -> bool;

        auto pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType;
        auto pattern_graph_row(int g, int p) const -> const SVOBitset &;
        auto target_graph_row(int g, int t) const -> const SVOBitset &;

        auto forward_target_graph_row(int t) const -> const SVOBitset &;
        auto reverse_target_graph_row(int t) const -> const SVOBitset &;

        auto pattern_degree(int g, int p) const -> unsigned;
        auto target_degree(int g, int t) const -> unsigned;
        auto largest_target_degree() const -> unsigned;

        auto has_vertex_labels() const -> bool;
        auto has_edge_labels() const -> bool;
        auto directed() const -> bool;
        auto pattern_vertex_label(int p) const -> int;
        auto target_vertex_label(int p) const -> int;
        auto pattern_edge_label(int p, int q) const -> int;
        auto target_edge_label(int t, int u) const -> int;
        auto edge_label_compatibility(int p_lid, int t_lid) const -> bool;

        auto pattern_has_loop(int p) const -> bool;
        auto target_has_loop(int t) const -> bool;

        auto initialise_domains(std::vector<HomomorphismDomain> & domains) const -> bool;
};

#endif
