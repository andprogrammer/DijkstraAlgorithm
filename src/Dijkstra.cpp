//Print minimum distances from source vertex to all others
//Solution : Dijkstras shortest path algorithm (Greedy Algorithm)
#include "LibsAndClassDeclarations.h"


namespace Dijkstra
{
	using SHORTES_VERTEX_DISTANCE_TYPE = std::unordered_map<uint32_t, std::pair<uint32_t, bool>>;
	using GRAPH_TYPE = std::unordered_map<uint32_t, std::forward_list<std::pair<uint32_t, uint32_t>>>;
	auto MAX_UINT_32_T = std::numeric_limits<uint32_t>::max();

	void showGraph(const GRAPH_TYPE &graph)
	{
		std::cout << "\nGRAPH vertex - [adjacent : weight]\n" << std::endl;
		for (const auto &vertex : graph)
		{
			std::cout << vertex.first << " - ";
			for (const auto &adjacent : vertex.second)
			{
				std::cout << "[" << adjacent.first << " : " << adjacent.second << "] ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	void fillShortestVertexDistance(const GRAPH_TYPE &graph, SHORTES_VERTEX_DISTANCE_TYPE &shortestVertexDistances, uint32_t sourceVertex)
	{
		for (const auto &vertex : graph)
			shortestVertexDistances.insert({ vertex.first, {MAX_UINT_32_T, false} });	//std::initilizer_list faster solution in that case construct object - std::make_pair<uint32_t, bool>() slower - construct object, then move object

		shortestVertexDistances[sourceVertex].first = 0u;
	}

	uint32_t getMinimumDistance(const SHORTES_VERTEX_DISTANCE_TYPE &shortestVertexDistances)
	{
		uint32_t minimumDistance = MAX_UINT_32_T, minimum_vertex = 0u;
		for (const auto &vertex : shortestVertexDistances)
		{
			if (vertex.second.second == false && vertex.second.first <= minimumDistance)
				minimumDistance = vertex.second.first, minimum_vertex = vertex.first;
		}
		return minimum_vertex;
	}

	bool ifPathExist(const GRAPH_TYPE &graph, int32_t sourceVertex, int32_t destinationVertex)
	{
		for (const auto &vertex : graph)
		{
			if (vertex.first == sourceVertex)
			{
				for (const auto &adjacent : vertex.second)
				{
					if (adjacent.first == destinationVertex) return true;
				}
				return false;
			}
		}
		return false;
	}

	uint32_t edgeWeight(const GRAPH_TYPE &graph, uint32_t sourceVertex, uint32_t destinationVertex)
	{
		for (const auto &vertex : graph)
		{
			if (vertex.first == sourceVertex)
			{
				for (const auto &adjacent : vertex.second)
				{
					if (adjacent.first == destinationVertex) return adjacent.second;
				}
				return 0u;
			}
		}
		return 0u;
	}

	void updateDistanceOfAdjacentVerticesOfThePickedVertex(SHORTES_VERTEX_DISTANCE_TYPE &shortestVertexDistances, std::uint32_t minimumDistance, const GRAPH_TYPE &graph)
	{
		for (auto &vertex : shortestVertexDistances)
		{
			bool pathExist = ifPathExist(graph, minimumDistance, vertex.first);
			uint32_t edgeDistance = edgeWeight(graph, minimumDistance, vertex.first);

			if (!vertex.second.second
				&& pathExist
				&& shortestVertexDistances[minimumDistance].first != MAX_UINT_32_T
				&& shortestVertexDistances[minimumDistance].first + edgeDistance < vertex.second.first)
				vertex.second.first = shortestVertexDistances[minimumDistance].first + edgeDistance;
		}
	}

	void findShortestPathForAllVertices(const GRAPH_TYPE &graph, SHORTES_VERTEX_DISTANCE_TYPE &shortestVertexDistances)
	{
		size_t graphSize = graph.size();
		uint32_t visitedCounter = 0u;
		while (visitedCounter < graphSize - 1)
		{
			std::uint32_t minimumDistance = getMinimumDistance(shortestVertexDistances);
			shortestVertexDistances[minimumDistance].second = true;

			updateDistanceOfAdjacentVerticesOfThePickedVertex(shortestVertexDistances, minimumDistance, graph);
			++visitedCounter;
		}
	}

	void showAllVertexDistanceFromSourceVertex(const SHORTES_VERTEX_DISTANCE_TYPE &shortestVertexDistances, uint32_t sourceVertex)
	{
		std::cout << "Dijkstras shortest path algorithm (Greedy Algorithm)" << std::endl;
		std::cout << "Vertex : Distance from source " << sourceVertex << std::endl << std::endl;
		for (const auto &vertex : shortestVertexDistances)
			std::cout << vertex.first << " : " << vertex.second.first << std::endl;
	}

	void findShortestPath(const GRAPH_TYPE &graph, uint32_t sourceVertex)
	{
		SHORTES_VERTEX_DISTANCE_TYPE shortestVertexDistances;
		fillShortestVertexDistance(graph, shortestVertexDistances, sourceVertex);

		findShortestPathForAllVertices(graph, shortestVertexDistances);

		showAllVertexDistanceFromSourceVertex(shortestVertexDistances, sourceVertex);
	}
}

namespace DijkstraUnitTests
{
	void runUnitTests(const Dijkstra::GRAPH_TYPE &graph, uint32_t sourceVertex)
	{
		std::cout << "\n############ Unit tests ############" << std::endl;
		std::cout << std::boolalpha;

		Dijkstra::SHORTES_VERTEX_DISTANCE_TYPE shortestVertexDistances;
		Dijkstra::fillShortestVertexDistance(graph, shortestVertexDistances, sourceVertex);
		std::cout << (shortestVertexDistances.size() == 9) << std::endl;

		std::cout << (Dijkstra::getMinimumDistance(shortestVertexDistances) == sourceVertex) << std::endl;

		std::cout << (Dijkstra::ifPathExist(graph, 4u, 3u) == true) << std::endl;
		std::cout << (Dijkstra::ifPathExist(graph, 8u, 0u) == false) << std::endl;

		std::cout << (Dijkstra::edgeWeight(graph, 2u, 3u) == 4u) << std::endl;
		std::cout << (Dijkstra::edgeWeight(graph, 6u, 3u) == 0u) << std::endl;

		Dijkstra::updateDistanceOfAdjacentVerticesOfThePickedVertex(shortestVertexDistances, sourceVertex, graph);
		std::cout << (shortestVertexDistances[1u].first == 5u) << std::endl;
		std::cout << (shortestVertexDistances[2u].first == 7u) << std::endl;
		std::cout << (shortestVertexDistances[3u].first == Dijkstra::MAX_UINT_32_T) << std::endl;

		Dijkstra::findShortestPathForAllVertices(graph, shortestVertexDistances);
		std::cout << (shortestVertexDistances[1u].first == 5u) << std::endl;
		std::cout << (shortestVertexDistances[2u].first == 6u) << std::endl;
		std::cout << (shortestVertexDistances[3u].first == 10u) << std::endl;

		Dijkstra::findShortestPath(graph, sourceVertex);
		std::cout << (shortestVertexDistances[0u].first == 0u) << std::endl;
		std::cout << (shortestVertexDistances[1u].first == 5u) << std::endl;
		std::cout << (shortestVertexDistances[2u].first == 6u) << std::endl;
		std::cout << (shortestVertexDistances[3u].first == 10u) << std::endl;
		std::cout << (shortestVertexDistances[4u].first == 12u) << std::endl;
		std::cout << (shortestVertexDistances[5u].first == 8u) << std::endl;
		std::cout << (shortestVertexDistances[6u].first == 13u) << std::endl;
		std::cout << (shortestVertexDistances[7u].first == 9u) << std::endl;
		std::cout << (shortestVertexDistances[8u].first == 18u) << std::endl;
	}
}

int main()
{
	const Dijkstra::GRAPH_TYPE graph =
	{
		{ 0u, { { 1u, 5u }, { 2u, 7u } } },
		{ 1u, { { 0u, 5u }, { 2u, 1u }, { 4u, 7u } } },
		{ 2u, { { 0u, 7u }, { 1u, 1u }, { 3u, 4u }, { 5u, 2u } } },
		{ 3u, { { 2u, 4u }, { 4u, 5u }, { 5u, 6u } } },
		{ 4u, { { 1u, 7u }, { 3u, 5u }, { 6u, 6u }, { 7u, 4u } } },
		{ 5u, { { 2u, 2u }, { 3u, 6u }, { 6u, 5u }, { 7u, 1u } } },
		{ 6u, { { 4u, 6u }, { 5u, 5u }, { 8u, 5u } } },
		{ 7u, { { 4u, 4u }, { 5u, 1u }, { 8u, 9u } } },
		{ 8u, { { 6u, 5u }, { 7u, 9u } } }
	};

	Dijkstra::showGraph(graph);

	auto sourceVertex = 0u;
	Dijkstra::findShortestPath(graph, sourceVertex);

	// Simple Unit Tests
	DijkstraUnitTests::runUnitTests(graph, sourceVertex);

	return 0;
}

//clear && mkdir build && cd build && cmake .. && make && ./Dijkstra
