/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// #include "base_homomorphism_model.hh"
#include "homomorphism_model.hh"
#include "homomorphism_traits.hh"
#include "configuration.hh"

#include <functional>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::greater;
using std::list;
using std::map;
using std::max;
using std::optional;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::to_string;
using std::vector;
using std::multiset;

// // TODO: Probably have an issue in here w inheritance.
// struct HomomorphismModel::Imp
// {
//     const HomomorphismParams & params;
//
//     vector<PatternAdjacencyBitsType> pattern_adjacencies_bits;
//     vector<SVOBitset> pattern_graph_rows;
//     vector<SVOBitset> target_graph_rows, forward_target_graph_rows, reverse_target_graph_rows;
//
//     vector<vector<int> > patterns_degrees, targets_degrees;
//     vector<vector<bool> > _edge_label_compatibility;
//     int largest_target_degree = 0;
//     // TODO: Change has_less_thans variable name.
//     bool has_less_thans = false, directed = false;
//
//     vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;
//     // TODO: Discuss whether this is better as a vector or a map.
//     // The map below is essentially a copy of InputGraph.edges. Should we just modify that directly?
//     // map<pair<int,int>, int> pattern_edge_labels, target_edge_labels;
//     vector<int> pattern_loops, target_loops;
//
//     vector<string> pattern_vertex_proof_names, target_vertex_proof_names;
//
//     Imp(const HomomorphismParams & p) :
//         params(p)
//     {
//     }
// };

HomomorphismModel::HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    BaseHomomorphismModel(target, pattern, params)
    // _imp(new Imp(params)),
    // max_graphs(calculate_n_shape_graphs(params)),
    // pattern_size(pattern.size()),
    // target_size(target.size())
{
    // Map each unique edge_label to an integer.
    map<multiset<string>, int> pattern_edge_labels_map; //, target_edge_labels_map;
    // Resize vector recording integers corresponding to each edge's label.
    // _imp->pattern_edge_labels.resize(pattern_size * pattern_size);
    _imp->target_edge_labels.resize(target_size * target_size);

    // Fill edge_labels_map labels -> int and edge_labels with labels.
    _record_edge_labels(pattern_edge_labels_map, pattern, _imp->pattern_edge_labels);

    // Record only the edges whose labels occur in the pattern graph.
    // TODO: Implement.
    _record_edge_labels(target_edge_labels_map, target, _imp->target_edge_labels);

    // // TODO: Find better way to get indices.
    // Form edge compatibility matrix.
    _imp->_edge_label_compatibility.resize(pattern_edge_labels_map.size(), vector<bool>(target_edge_labels_map.size()));
    for (const auto& [labels1, label1_id] : pattern_edge_labels_map) {
        for (const auto& [labels2, label2_id] : target_edge_labels_map) {
            // TODO: Deal with the induced case.
            _imp->_edge_label_compatibility[label1_id][label2_id] = check_edge_label_compatibility(labels1, labels2);
        }
    }

    // recode target to a bit graph, and take out loops
    _imp->target_graph_rows.resize(target_size * max_graphs, SVOBitset{ target_size, 0 });
    _imp->target_loops.resize(target_size);
    for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
        if (e->first.first == e->first.second)
            _imp->target_loops[e->first.first] = 1;
        else
            _imp->target_graph_rows[e->first.first * max_graphs + 0].set(e->first.second);
    }

    // if directed, do both directions
    if (pattern.directed()) {
        _imp->forward_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        _imp->reverse_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        // Record non-loopy edges.
        for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
            if (e->first.first != e->first.second) {
                _imp->forward_target_graph_rows[e->first.first].set(e->first.second);
                _imp->reverse_target_graph_rows[e->first.second].set(e->first.first);
            }
        }
    }

}

/** Check that at least as many of each label occur for the target as pattern edge.
TODO: Change to exact count for induced case.
*/
auto HomomorphismModel::check_edge_label_compatibility(const multiset<string>& labels1, const multiset<string>& labels2) const -> bool
{
    // Map labels to the number of occurences for each multiset.
    map<string, int> label_counts1 = _multiset_item_counts(labels1);
    map<string, int> label_counts2 = _multiset_item_counts(labels2);

    // Check compatibility for each label.
    for (const auto & [label, count] : label_counts1) {
        if (label_counts2.find(label) == label_counts2.end() || count > label_counts2[label])
            return false;
    }
    return true;
}

auto HomomorphismModel::edge_label_compatibility(int p_lid, int t_lid) const -> bool
{
    return _imp->_edge_label_compatibility[p_lid][t_lid];

}




        // // Print labels and compatibility matrix.
        // std::cout << "Pattern edge labels: " << std::endl;
        // for (const auto& [labels1, label1_id] : pattern_edge_labels_map) {
        //     std::cout << "Label set id: " << label1_id << ", label set: " << std::endl;
        //     for (const auto& label : labels1)
        //         std::cout << label << " ";
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
        // std::cout << "Target edge labels: " << std::endl;
        // for (const auto& [labels2, label2_id] : target_edge_labels_map) {
        //     std::cout << "Label set id: " << label2_id << ", label set: " << std::endl;
        //     for (const auto& label : labels2)
        //         std::cout << label << " ";
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
        //
        // std::cout << "Label compatibility: " << std::endl;
        // for (const auto& [labels1, label1_id] : pattern_edge_labels_map) {
        //     for (const auto& [labels2, label2_id] : target_edge_labels_map){
        //         std::cout << _imp->edge_label_compatibility[label1_id][label2_id] << " ";
        //     }
        //     std::cout << std::endl;
        // }
