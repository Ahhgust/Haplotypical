#ifndef GRAPHLINE_TYPES_H
#define GRAPHLINE_TYPES_H 1

#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/**A node in a sequence graph.*/
typedef struct{
	/**The number of backwards links.*/
	uintptr_t numPrev;
	/**The characters for those links.*/
	char* linkChars;
	/**The nodes for those links.*/
	uintptr_t* linkNodes;
} GraphlineNode;

/**A node in a sequence graph.*/
typedef struct{
	/**The number of backwards links.*/
	uintptr_t numPrev;
	/**The offset into the character array for those links.*/
	uintptr_t charOffset;
	/**The offset into the node array for those links.*/
	uintptr_t nodeOffset;
} GraphlineNodeStore;

/**A sequence graph.*/
class GraphlineGraph{
public:
	
	/**Set up a graph with a single start node.*/
	GraphlineGraph();
	/**Tear down.*/
	~GraphlineGraph();
	
	/**
	 * Returns the number of nodes in this graph.
	 * @param The number of nodes in this graph.
	 */
	uintptr_t size();
	
	/**
	 * Get a node from this graph.
	 * @param nodeInd The index to get.
	 * @param toFill The place to put the data.
	 */
	void getNode(uintptr_t nodeInd, GraphlineNode* toFill);
	
	/**
	 * Add a node to this graph.
	 * @param numLink THe number of backwards links.
	 * @param linkChars The characters for those links.
	 * @param linkNodes The nodes for those links.
	 */
	void addNode(int numLink, const char* linkChars, uintptr_t* linkNodes);
	
	/**
	 * Return the topological sorting of this graph (all nodes with no inputs first)..
	 * @param toFill The place to put the node indices.
	 */
	void topologicalSort(std::vector<uintptr_t>* toFill);
	
	/**All the nodes.*/
	std::vector<GraphlineNodeStore> allNodes;
	/**The number of links baked in.*/
	uintptr_t numBake;
	/**The capacity of said links.*/
	uintptr_t capBake;
	/**The characters of the links.*/
	char* linkCharArr;
	/**The nodes of the links.*/
	uintptr_t* linkNodeArr;
	
	/**Whether the topological sort is still valid.*/
	bool topoValid;
	/**The topological sort.*/
	std::vector<uintptr_t> lastSort;
	/**The nodes that lead nowhere.*/
	std::vector<uintptr_t> termNodes;
};

/**Stores the result of an alignment.*/
class AlignmentResult{
public:
	/**
	 * Stores the result of an alignment.
	 * @param alignTable The alignment table.
	 * @param numNodes The number of nodes in the original graph.
	 * @param seqLength The length of the sequence being aligned.
	 */
	AlignmentResult(int** alignTable, std::string* origS, GraphlineGraph* forGrap);
	/**Clean up.*/
	~AlignmentResult();
	/**The number of nodes in the original graph.*/
	int numNode;
	/**The alignment table.*/
	int** alnTab;
	/**The terminal nodes of the graph.*/
	std::vector<uintptr_t> endpoints;
	/**The original string.*/
	std::string origString;
};




#endif
