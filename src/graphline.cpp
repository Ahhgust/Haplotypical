
#include "Haplotypical_types.h"

/**
 * This will calculate the edit distance table.
 * @param stringLen The length of the string.
 * @param forString The string.
 * @param forGraph The graph to align to.
 * @return The table. First index is point in string (0 is before the first character), second index is the node index.
 */
int** calculateStringGraphAlignmentTable(int stringLen, const char* forString, GraphlineGraph* forGraph);

/**
 * Traverse an alignment.
 * @param stringLen The length of the strings.
 * @param forString The string.
 * @param forGraph The graph to work over.
 * @param alnTab The alignment table.
 * @param fillS The string indices.
 * @param fillG The graph node indices.
 * @return The number of filled indices.
 */
int traverseGraphAlignTable(int stringLen, const char* forString, GraphlineGraph* forGraph, int** alnTab, int* fillS, int* fillG);

//To put in its own source

#include <map>
#include <iostream>
#include <algorithm>

GraphlineGraph::GraphlineGraph(){
	topoValid = false;
	numBake = 0;
	capBake = 64;
	linkCharArr = (char*)malloc(capBake);
	linkNodeArr = (uintptr_t*)malloc(capBake*sizeof(uintptr_t));
	//add the start node
	addNode(0, 0, 0);
}

GraphlineGraph::~GraphlineGraph(){
	free(linkCharArr);
	free(linkNodeArr);
}

uintptr_t GraphlineGraph::size(){
	return allNodes.size();
}

void GraphlineGraph::getNode(uintptr_t nodeInd, GraphlineNode* toFill){
	GraphlineNodeStore tmpFill = allNodes[nodeInd];
	toFill->numPrev = tmpFill.numPrev;
	toFill->linkChars = linkCharArr + tmpFill.charOffset;
	toFill->linkNodes = linkNodeArr + tmpFill.nodeOffset;
}

void GraphlineGraph::addNode(int numLink, const char* linkChars, uintptr_t* linkNodes){
	uintptr_t needSize = numLink + numBake;
	if(needSize > capBake){
		while(needSize > capBake){
			capBake = 2*capBake;
		}
		char* nextCharArr = (char*)malloc(capBake);
		uintptr_t* nextNodeArr = (uintptr_t*)malloc(capBake*sizeof(uintptr_t));
		memcpy(nextCharArr, linkCharArr, numBake);
		memcpy(nextNodeArr, linkNodeArr, numBake*sizeof(uintptr_t));
		free(linkCharArr);
		free(linkNodeArr);
		linkCharArr = nextCharArr;
		linkNodeArr = nextNodeArr;
	}
	GraphlineNodeStore tmpFill;
		tmpFill.numPrev = numLink;
		tmpFill.charOffset = numBake;
		tmpFill.nodeOffset = numBake;
	allNodes.push_back(tmpFill);
	memcpy(linkCharArr + numBake, linkChars, numLink);
	memcpy(linkNodeArr + numBake, linkNodes, numLink*sizeof(uintptr_t));
	numBake += numLink;
	topoValid = false;
}

void GraphlineGraph::topologicalSort(std::vector<uintptr_t>* toFill){
	if(!topoValid){
		lastSort.clear();
		termNodes.clear();
		//swapping vectors
		std::vector<uintptr_t> curTakeNode;
		std::vector<uintptr_t> nextTakeNode;
		//build a working set of info (number of incoming, and where outgoing goes)
		std::map<uintptr_t, std::vector<uintptr_t> > outgoingMap;
		std::vector<uintptr_t> numIncome;
		for(uintptr_t i = 0; i<allNodes.size(); i++){
			GraphlineNodeStore tmpFill = allNodes[i];
			numIncome.push_back(tmpFill.numPrev);
			for(uintptr_t j = 0; j<tmpFill.numPrev; j++){
				uintptr_t fromNode = linkNodeArr[tmpFill.nodeOffset + j];
				outgoingMap[fromNode].push_back(i);
			}
			if(tmpFill.numPrev == 0){
				curTakeNode.push_back(i);
			}
		}
		//remove zero nodes until the graph is empty
		while(curTakeNode.size()){
			nextTakeNode.clear();
			for(uintptr_t i = 0; i<curTakeNode.size(); i++){
				uintptr_t curTake = curTakeNode[i];
				lastSort.push_back(curTake);
				std::vector<uintptr_t>* nodeUnlink = &(outgoingMap[curTake]);
				for(uintptr_t j = 0; j<nodeUnlink->size(); j++){
					uintptr_t curUL = (*nodeUnlink)[j];
					numIncome[curUL] = numIncome[curUL] - 1;
					if(numIncome[curUL] == 0){
						nextTakeNode.push_back(curUL);
					}
				}
			}
			curTakeNode.clear();
			curTakeNode.insert(curTakeNode.begin(), nextTakeNode.begin(), nextTakeNode.end());
		}
		//note which nodes have no outgoing
		for(uintptr_t i = 0; i<allNodes.size(); i++){
			if(outgoingMap[i].size() == 0){
				termNodes.push_back(i);
			}
		}
	}
	toFill->clear();
	toFill->insert(toFill->end(), lastSort.begin(), lastSort.end());
}

int** calculateStringGraphAlignmentTable(int stringLen, const char* forString, GraphlineGraph* forGraph){
	//set up the table pointers
		int** fullTable = (int**)malloc((stringLen+1)*sizeof(int*) + (stringLen+1)*(forGraph->size())*sizeof(int));
		int* curFoc = (int*)(fullTable + (stringLen+1));
		for(int i = 0; i<=stringLen; i++){
			fullTable[i] = curFoc;
			curFoc += forGraph->size();
		}
	//note the topological sorting
		std::vector<uintptr_t> nodeSort;
		forGraph->topologicalSort(&nodeSort);
	//run along the string
		GraphlineNode focNode;
		for(int i = 0; i<=stringLen; i++){
			//int* curLine = fullTable[i];
			//run down the nodes
			for(uintptr_t j = 0; j<nodeSort.size(); j++){
				uintptr_t curN = nodeSort[j];
				forGraph->getNode(curN, &focNode);
				//find the minimum distance
				bool minSkip = true;
				int curMin = 0;
				int testMin;
				#define RUN_NEW_SCORE(newSEq) \
					testMin = newSEq;\
					if(minSkip){ curMin = testMin; minSkip = false; }\
					else if(testMin < curMin){ curMin = testMin; }
				if(i){
					RUN_NEW_SCORE(fullTable[i-1][curN] + 1);
					for(uintptr_t k = 0; k<focNode.numPrev; k++){
						RUN_NEW_SCORE(fullTable[i-1][focNode.linkNodes[k]] + (focNode.linkChars[k] == forString[i-1] ? 0 : 1))
					}
				}
				for(uintptr_t k = 0; k<focNode.numPrev; k++){
					RUN_NEW_SCORE(fullTable[i][focNode.linkNodes[k]] + 1)
				}
				fullTable[i][curN] = curMin;
			}
		}
	return fullTable;
}

int traverseGraphAlignTable(int stringLen, const char* forString, GraphlineGraph* forGraph, int** alnTab, int* fillS, int* fillG){
	GraphlineNode focNode;
	//start at the terminal node with the lowest score
	uintptr_t curS = stringLen;
	uintptr_t curN = forGraph->termNodes[0];
	int curNScr = alnTab[curS][curN];
	for(unsigned i = 1; i<forGraph->termNodes.size(); i++){
		uintptr_t testN = forGraph->termNodes[i];
		int testScr = alnTab[curS][testN];
		if(testScr < curNScr){
			curN = testN;
			curNScr = testScr;
		}
	}
	*fillS = curS; *fillG = curN;
	int numAdd = 1;
	while(curN || curS){
		forGraph->getNode(curN, &focNode);
		int curMin = alnTab[curS][curN];
		int testMin;
		int winN = -1; int winS = -1;
		#define TEST_NEW_SCORE(newSEq, testN, testS) \
			testMin = newSEq;\
			if(testMin == curMin){ winN = testN; winS = testS; }
		if(curS){
			TEST_NEW_SCORE(alnTab[curS-1][curN] + 1, curN, curS - 1)
			for(uintptr_t k = 0; k<focNode.numPrev; k++){
				TEST_NEW_SCORE(alnTab[curS-1][focNode.linkNodes[k]] + (focNode.linkChars[k] == forString[curS-1] ? 0 : 1), focNode.linkNodes[k], curS - 1)
			}
		}
		for(uintptr_t k = 0; k<focNode.numPrev; k++){
			TEST_NEW_SCORE(alnTab[curS][focNode.linkNodes[k]] + 1, focNode.linkNodes[k], curS)
		}
		curN = winN; curS = winS;
		fillS[numAdd] = curS;
		fillG[numAdd] = curN;
		numAdd++;
	}
	std::reverse(fillS, fillS + numAdd);
	std::reverse(fillG, fillG + numAdd);
	return numAdd;
}


//RCpp stuff

#include <Rcpp.h>
#include <iostream>

/**A node to use while building the graph.*/
class BuildGraphNode{
public:
	/**The sequence of this node.*/
	std::string nodeSeq;
	/**The nodes this points at.*/
	std::vector<uintptr_t> nextNodes;
};

/**The info on the event.*/
class EventInfo{
public:
	/**The allele.*/
	std::string allele;
	/**The position.*/
	int position;
	/**The type.*/
	int evtType;
};

/**
 * Calculates forward links, taking epsilon jumps into account.
 * @param allNodes All the nodes.
 * @param fromNI The index to work over.
 * @param allFound The link index vectors to fill in.
 */
void getAllForwardWithEpsilon(std::vector<BuildGraphNode>* allNodes, uintptr_t fromNI, std::vector< std::vector<uintptr_t> >* allFound){
	//quick stop if already calculated
		std::vector<uintptr_t>* curFocN = &((*allFound)[fromNI]);
		if(curFocN->size()){
			return;
		}
		BuildGraphNode* curFoc = &((*allNodes)[fromNI]);
	//get the full run for all the nexts
		for(unsigned i = 0; i<curFoc->nextNodes.size(); i++){
			getAllForwardWithEpsilon(allNodes, curFoc->nextNodes[i], allFound);
		}
	//add each to the current (or full if is an epsilon)
		for(unsigned i = 0; i<curFoc->nextNodes.size(); i++){
			uintptr_t focI = curFoc->nextNodes[i];
			if((*allNodes)[focI].nodeSeq.size()){
				curFocN->push_back(curFoc->nextNodes[i]);
			}
			else{
				std::vector<uintptr_t>* nextFoc = &((*allFound)[focI]);
				curFocN->insert(curFocN->end(), nextFoc->begin(), nextFoc->end());
			}
		}
}

#define EVENT_TYPE_BASE 1
#define EVENT_TYPE_INSERT 2
#define EVENT_TYPE_DELETE 3



/**
 * Compare two events.
 * @param coma The first thing to compare.
 * @param comb The second thing to compare.
 * @return Whether coma is before comb.
 */
bool eventCompare(const EventInfo& coma, const EventInfo& comb){
	return coma.position < comb.position;
}



//' Given a description of the variants in a mixture, will generate a graph for use with later analysis.
//' 
//' @param refRS The reference sequence.
//' @param diffAllele The variant alleles for each difference.
//' @param position The positions in the reference: one based. Multiple events can be specified at a location.
//' @param etype The type of each event. Zero for replacement, one for insert, two for delete.
//' @return An opaque object for use with other methods in this package.
// [[Rcpp::export]]
Rcpp::XPtr<GraphlineGraph> makeSequenceGraph(Rcpp::String refRS, Rcpp::StringVector diffAllele, Rcpp::IntegerVector position, Rcpp::IntegerVector etype){
	std::string ref = refRS;
	const char* refStr = ref.c_str();
	//idiot checks
	if(diffAllele.size() != position.size()){
		throw Rcpp::exception("Input vector sizes do not match.");
	}
	if(diffAllele.size() != etype.size()){
		throw Rcpp::exception("Input vector sizes do not match.");
	}
	for(int i = 0; i<diffAllele.size(); i++){
		if(position[i] < 0){
			throw Rcpp::exception("Cannot have an event at a negative index.");
		}
		else if((position[i] == 0) && (etype[i] != EVENT_TYPE_INSERT)){
			throw Rcpp::exception("Can only have insertion events before the start of the string.");
		}
	}
	//sort the variants by location
		std::vector<EventInfo> allVarInf;
		for(unsigned i = 0; i<diffAllele.size(); i++){
			EventInfo curAdd;
			curAdd.allele = diffAllele[i];
			curAdd.position = position[i];
			curAdd.evtType = etype[i];
			allVarInf.push_back(curAdd);
		}
		std::sort(allVarInf.begin(), allVarInf.end(), eventCompare);
	//build up a graph with direct links
		int curRI = ref.size() - 1;
		std::vector<BuildGraphNode> allNodes;
		std::vector<uintptr_t> curHorizon;
		int curEI = allVarInf.size() - 1;
		while(curEI >= 0){
			//find all events at the same position
				int endEI = curEI;
				int endPos = allVarInf[endEI].position;
				if(endPos == 0){
					break;
				}
				curEI--;
				while(curEI >= 0){
					if(allVarInf[curEI].position != endPos){
						break;
					}
					curEI--;
				}
				endPos--; //back to zero based
			//add the reference up to this location
				if(curRI > endPos){
					BuildGraphNode seqNode;
					seqNode.nodeSeq = std::string(refStr + endPos + 1, refStr + curRI + 1);
					seqNode.nextNodes = curHorizon;
					allNodes.push_back(seqNode);
					curHorizon.clear();
					curHorizon.push_back(allNodes.size()-1);
				}
			//add any insertions
				bool anyInsert = false;
				for(int i = curEI + 1; i <= endEI; i++){
					if(allVarInf[i].evtType == EVENT_TYPE_INSERT){
						anyInsert = true;
						break;
					}
				}
				if(anyInsert){
					std::vector<uintptr_t> nextHorizon;
					BuildGraphNode skipHor;
						skipHor.nextNodes = curHorizon;
						allNodes.push_back(skipHor);
						nextHorizon.push_back(allNodes.size()-1);
					for(int i = curEI + 1; i <= endEI; i++){
						if(allVarInf[i].evtType != EVENT_TYPE_INSERT){
							continue;
						}
						BuildGraphNode insertHor;
						insertHor.nodeSeq = allVarInf[i].allele;
						insertHor.nextNodes = curHorizon;
						allNodes.push_back(insertHor);
						nextHorizon.push_back(allNodes.size()-1);
					}
					curHorizon = nextHorizon;
				}
			//add all variants (including deletions)
				std::vector<uintptr_t> nextHorizon;
				BuildGraphNode refHor;
					refHor.nodeSeq.push_back(refStr[endPos]);
					refHor.nextNodes = curHorizon;
					allNodes.push_back(refHor);
					nextHorizon.push_back(allNodes.size()-1);
				for(int i = curEI + 1; i <= endEI; i++){
					if(allVarInf[i].evtType == EVENT_TYPE_INSERT){
						continue;
					}
					BuildGraphNode insertHor;
					if(allVarInf[i].evtType == EVENT_TYPE_BASE){
						insertHor.nodeSeq = allVarInf[i].allele;
					}
					insertHor.nextNodes = curHorizon;
					allNodes.push_back(insertHor);
					nextHorizon.push_back(allNodes.size()-1);
				}
				curHorizon = nextHorizon;
			curRI = endPos - 1;
		}
	//finish off the reference
		if(curRI >= 0){
			BuildGraphNode seqNode;
			seqNode.nodeSeq = std::string(refStr, refStr + curRI + 1);
			seqNode.nextNodes = curHorizon;
			allNodes.push_back(seqNode);
			curHorizon.clear();
			curHorizon.push_back(allNodes.size()-1);
		}
	//add any remaining insertions, if any
		if(curEI >= 0){
			std::vector<uintptr_t> nextHorizon;
			BuildGraphNode skipHor;
				skipHor.nextNodes = curHorizon;
				allNodes.push_back(skipHor);
				nextHorizon.push_back(allNodes.size()-1);
			for(int i = 0; i <= curEI; i++){
				BuildGraphNode insertHor;
				insertHor.nodeSeq = allVarInf[i].allele;
				insertHor.nextNodes = curHorizon;
				allNodes.push_back(insertHor);
				nextHorizon.push_back(allNodes.size()-1);
			}
			curHorizon = nextHorizon;
		}
	//add a starting node
		uintptr_t startNI = allNodes.size();
		BuildGraphNode startNode;
			startNode.nextNodes = curHorizon;
			allNodes.push_back(startNode);
	//note the full extent of the forward links (of the simple nodes)
		std::vector< std::vector<uintptr_t> > allForward;
		allForward.resize(allNodes.size());
		for(unsigned i = 0; i<allNodes.size(); i++){
			getAllForwardWithEpsilon(&allNodes, i, &allForward);
		}
	//turn them into backward links
		std::vector< std::vector<uintptr_t> > allBackward;
		allBackward.resize(allNodes.size());
		for(unsigned i = 0; i<allNodes.size(); i++){
			if(allNodes[i].nodeSeq.size() == 0){
				continue;
			}
			std::vector<uintptr_t>* curFL = &(allForward[i]);
			for(unsigned j = 0; j<curFL->size(); j++){
				uintptr_t fromI = (*curFL)[j];
				std::vector<uintptr_t>* addBack = &(allBackward[fromI]);
				addBack->push_back(i);
			}
		}
	//note which items need a back link to the start node
		std::set<uintptr_t> needBackStart;
		{
			std::vector<uintptr_t>* curFL = &(allForward[startNI]);
			needBackStart.insert(curFL->begin(), curFL->end());
		}
	//turn the direct graph into the full graph (start at the horizon and work out)
		GraphlineGraph* newGraph = new GraphlineGraph();
		std::map<uintptr_t,uintptr_t> nodeBackMap; //from BuildGraphNode index to newGraph index
		std::set<uintptr_t> needHand;
			//start by looking at the targets of the first node
			needHand.insert(allForward[startNI].begin(), allForward[startNI].end());
		//some storage space
			std::vector<uintptr_t> backLinkRT;
			std::string backLinkCs;
		//loop until all nodes handled
		while(needHand.size()){
			std::set<uintptr_t> nextHand;
			for(std::set<uintptr_t>::iterator curLook = needHand.begin(); curLook != needHand.end(); curLook++){
				uintptr_t curLI = *curLook;
				//skip if already handled
					if(nodeBackMap.find(curLI) != nodeBackMap.end()){
						continue;
					}
				//make sure all back links are handled
					bool haveAllBL = true;
					std::vector<uintptr_t>* curBack = &(allBackward[curLI]);
					for(unsigned i = 0; i<curBack->size(); i++){
						if(nodeBackMap.find((*curBack)[i]) == nodeBackMap.end()){
							haveAllBL = false;
							break;
						}
					}
					if(!haveAllBL){
						continue;
					}
				//clean up and build the links
					backLinkRT.clear();
					backLinkCs.clear();
					std::string* curNodS = &(allNodes[curLI].nodeSeq);
					for(unsigned i = 0; i<curBack->size(); i++){
						backLinkRT.push_back(nodeBackMap[(*curBack)[i]]);
						backLinkCs.push_back((*curNodS)[0]);
					}
					if(needBackStart.count(curLI)){
						backLinkRT.push_back(0);
						backLinkCs.push_back((*curNodS)[0]);
					}
				//add the first node
					uintptr_t* backLinkRTP = (backLinkRT.size() ? &(backLinkRT[0]) : 0);
					const char* backLinkCsP = (backLinkRT.size() ? &(backLinkCs[0]) : 0);
					newGraph->addNode(backLinkRT.size(), backLinkCsP, backLinkRTP);
					backLinkRT.clear();
					backLinkCs.clear();
				//add the rest of the nodes
					for(unsigned i = 1; i<curNodS->size(); i++){
						backLinkRT.push_back(newGraph->size()-1);
						backLinkCs.push_back((*curNodS)[i]);
						newGraph->addNode(backLinkRT.size(), &(backLinkCs[0]), &(backLinkRT[0]));
						backLinkRT.clear();
						backLinkCs.clear();
					}
				//note where it goes back to
					nodeBackMap[curLI] = (newGraph->size()-1);
				//add all the forward to the next handle
					std::vector<uintptr_t>* curFor = &(allForward[curLI]);
					nextHand.insert(curFor->begin(), curFor->end());
			}
			needHand = nextHand;
		}
	Rcpp::XPtr<GraphlineGraph> toRet(newGraph);
	return toRet;
}

AlignmentResult::AlignmentResult(int** alignTable, std::string* origS, GraphlineGraph* forGrap){
	numNode = forGrap->size();
	alnTab = alignTable;
	endpoints = forGrap->termNodes;
	origString = *origS;
}
AlignmentResult::~AlignmentResult(){
	free(alnTab);
}

//' Get the alignment table from a sequence to a graph.
//' 
//' @param sgrap The sequence graph in question.
//' @param query The sequence to align.
//' @return The alignment table.
// [[Rcpp::export]]
Rcpp::XPtr<AlignmentResult> alignSequenceGraph(Rcpp::XPtr<GraphlineGraph> sgrap, Rcpp::String query){
	std::string queryS = query;
	const char* quec = query.get_cstring();
	int** alnTab = calculateStringGraphAlignmentTable(strlen(quec), quec, sgrap.get());
	AlignmentResult* toRetV = new AlignmentResult(alnTab, &queryS, sgrap.get());
	Rcpp::XPtr<AlignmentResult> toRet(toRetV);
	return toRet;
}

//' Get the edit distance from a sequence to a graph.
//' 
//' @param alnTable The result of the alignment.
//' @return The edit distance between the graph and the query.
// [[Rcpp::export]]
int getSequenceGraphEditDistance(Rcpp::XPtr<AlignmentResult> alnTable){
	AlignmentResult* alnR = alnTable.get();
	uintptr_t curS = alnR->origString.size();
	int minV = alnR->alnTab[curS][alnR->endpoints[0]];
	for(unsigned i = 1; i<alnR->endpoints.size(); i++){
		int curV = alnR->alnTab[curS][alnR->endpoints[i]];
		if(curV < minV){
			minV = curV;
		}
	}
	return minV;
}

//' Note which nodes of the graph were not visited in an alignment.
//' 
//' @param sgrap The original graph.
//' @param alnTable The result of the alignment.
//' @return The nodes of the graph which were not visited.
// [[Rcpp::export]]
Rcpp::IntegerVector getAlignmentMissedNodes(Rcpp::XPtr<GraphlineGraph> sgrap, Rcpp::XPtr<AlignmentResult> alnTable){
	AlignmentResult* alnR = alnTable.get();
	GraphlineGraph* forGraph = sgrap.get();
	std::set<uintptr_t> visitedNodes;
	int* sinds = (int*)malloc((alnR->origString.size() + forGraph->size() + 2)*sizeof(int));
	int* ninds = (int*)malloc((alnR->origString.size() + forGraph->size() + 2)*sizeof(int));
	int numNV = traverseGraphAlignTable(alnR->origString.size(), alnR->origString.c_str(), forGraph, alnR->alnTab, sinds, ninds);
	for(int i = 0; i<numNV; i++){
		visitedNodes.insert(ninds[i]);
	}
	free(sinds);
	free(ninds);
	Rcpp::IntegerVector toRet;
	for(unsigned i = 0; i<forGraph->size(); i++){
		if(visitedNodes.find(i) == visitedNodes.end()){
			toRet.push_back(i);
		}
	}
	return toRet;
}

//' Prepare a string representation of the graph.
//' 
//' @param sgrap The sequence graph in question.
//' @return A string representation of the graph.
// [[Rcpp::export]]
std::string graphToString(Rcpp::XPtr<GraphlineGraph> sgrap){
	std::string toRet;
	char buffer[7 + 8*sizeof(int)+8];
	GraphlineGraph* forGraph = sgrap.get();
	for(unsigned i = 0; i<forGraph->size(); i++){
		GraphlineNode curNod;
		forGraph->getNode(i, &curNod);
		if(curNod.numPrev == 0){
			sprintf(buffer, "S %d\n", (int)i);
			toRet.append(buffer);
		}
		for(unsigned j = 0; j<curNod.numPrev; j++){
			char linkC = curNod.linkChars[j];
			uintptr_t preI = curNod.linkNodes[j];
			sprintf(buffer, "%c %d - %d\n", linkC, (int)preI, (int)i);
			toRet.append(buffer);
		}
	}
	return toRet;
}

