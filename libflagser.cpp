#include "flagser.hpp"

#define main no_dont_name_this_main
#include "flagser/src/flagser-count.cpp"
#undef main

#include <memory>

extern "C" size_t* flagser_count_unweighted(size_t nvertices, size_t nedges, const vertex_index_t edges[][2], size_t *res_size) {
	std::vector<value_t> vertex_filtration(nvertices, value_t(0));
	filtered_directed_graph_t graph(vertex_filtration, true);
	for (size_t edgei = 0; edgei < nedges; edgei ++){
		graph.add_filtered_edge(edges[edgei][0], edges[edgei][1], 0);
	}
	flagser_parameters params;
	params.nb_threads = 1;
	auto cout_buff = std::cout.rdbuf();
	std::cout.rdbuf(nullptr);
	auto res_cpp = count_cells(graph, params);
	std::cout.rdbuf(cout_buff);
	*res_size = res_cpp.size();
	void *res = malloc(sizeof(size_t) * res_cpp.size());
	memcpy(res, res_cpp.data(), res_cpp.size() * sizeof(size_t));
	return (size_t*)res;
}
