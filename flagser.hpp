#pragma once

#define MANY_VERTICES

#include "flagser/include/directed_graph.h"
#include "flagser/include/parameters.h"

#include <vector>

std::vector<size_t> count_cells(filtered_directed_graph_t& graph, const flagser_parameters& params);

extern "C" size_t* flagser_count_unweighted(size_t nvertices, size_t nedges, const vertex_index_t edges[][2], size_t *res_size);
