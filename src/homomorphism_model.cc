/* vim: set sw=4 sts=4 et foldmethod=syntax :

Extends BaseHomomorphismModel to perform precomputations on target edges.
Should be used when computationally feasible.
*/

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

HomomorphismModel::HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    BaseHomomorphismModel(target, pattern, params)
{
    // Resize vector recording integers corresponding to each edge's label.
    _imp->target_edge_labels.resize(target_size * target_size);

    // Fill edge_labels_map labels -> int and edge_labels with labels.
    _record_edge_labels(_imp->target_edge_labels_map, target, _imp->target_edge_labels);

    // // TODO: Find better way to get indices.
    // Form edge compatibility matrix.
    edge_label_compatibility.resize(_imp->pattern_edge_labels_map.size(), vector<bool>(_imp->target_edge_labels_map.size()));
    for (const auto& [labels1, label1_id] : _imp->pattern_edge_labels_map) {
        for (const auto& [labels2, label2_id] : _imp->target_edge_labels_map) {
            // TODO: Deal with the induced case.
            edge_label_compatibility[label1_id][label2_id] = check_edge_label_compatibility(labels1, labels2);
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

auto HomomorphismModel::check_edge_label_compatibility(int p_lid, int t_lid) const -> bool
{
    return edge_label_compatibility[p_lid][t_lid];
}

auto HomomorphismModel::target_edge_label(int t, int u) const -> int
{
    return _imp->target_edge_labels[t * target_size + u];
}
