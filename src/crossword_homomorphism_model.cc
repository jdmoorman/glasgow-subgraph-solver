
#include "crossword_homomorphism_model.hh"

#include <vector>
#include <string>
#include <set>
#include <stdexcept>

using std::vector;
using std::string;
using std::multiset;

CrosswordHomomorphismModel::CrosswordHomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    BaseHomomorphismModel(target, pattern, params)
{
    // Fill vector of vertex names.
    vertex_names.resize(target_size);
    for (unsigned i = 0 ; i < target_size ; ++i) {
        vertex_names[i] = target.vertex_name(i);
    }

    // Fill vector mapping ids to pattern edge label multisets.
    id_to_pattern_edge_labels.resize(_imp->pattern_edge_labels_map.size());
    for (auto const& [label, label_id] : _imp->pattern_edge_labels_map){
        id_to_pattern_edge_labels[label_id] = label;
    }
}

auto CrosswordHomomorphismModel::check_edge_label_compatibility(const int t_v1, const int t_v2, const int p_lid) const -> bool{

    multiset<string> p_label = id_to_pattern_edge_labels[p_lid];
    // Ensure there is a single label for the pattern edge.
    if (p_label.size() != 1) {
        throw std::invalid_argument( "More than one label for crossword pattern." );
    }

    // Extract intersection position.
    int i, j;
    for (multiset<string>::iterator it=p_label.begin(); it!=p_label.end(); ++it) {
        string label = *it;
        auto delim_pos = label.find('_');
        i = stoi(label.substr(0, delim_pos));
        j = stoi(label.substr(delim_pos + 1));
    }

    // Indicate whether the indicated intersection occurs in the vertices.
    return vertex_names[t_v1][i] == vertex_names[t_v2][j];
}
