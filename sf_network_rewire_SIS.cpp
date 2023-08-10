#include<iostream>
#include <algorithm>
#include <string>
#include <map>
#include <array>
#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<time.h>
#include<ctime>
#include<random>
#include<thread>
#include<sstream>
#include <list>
#include <iterator>
#include <vector>

//----
/*
Code to generate a Scale-Free Network following the Barabasi-Albert Algorithm. Also, it is possible to use a rewiring method
to reconnect nodes and study the propagation of a disease using the SIS epidemiological method.
*/
//----

using namespace std;

#define TAU 6.283185307


class Edge;
class Node;

class Edge
{
	public:
		int destinationID;
		int n_edges = 0;

		Edge() {}
		Edge(int destID){
			destinationID = destID;
		}

		int getDestinationID(){
			return destinationID;
		}

};

class Node {

public:
	int node_id;
	int node_status; //0=S, 1=I, 2=R
	int recovery_time = 10;
	int dead_time;
	list <Edge> edgeList;

	Node(){
		node_id = 0;
		node_status = 0;
		recovery_time = 10;
		dead_time = 2;
	}

	void printEdgeList() {
    	cout << "[";
    	for (auto it = edgeList.begin(); it != edgeList.end(); it++) {
      		cout << it -> getDestinationID() << " --> ";
    	}
    	cout << "]";
    	cout << endl;
    }

    void printNumberEdges(){
    	//cout << "Number of Edges for this node: " << edgeList.size() << endl;
    }


	void setID(int id){
		node_id = id;
	}

	int getID(){
		return node_id;
	}

	void setStatus(int status){
		node_status = status;
	}

	int getStatus(){
		return node_status;
	}

	list <Edge> getEdgeList() {
		return edgeList;
	}

	void addEdgeToEdgelist(int toNodeID)
    {
  		Edge e(toNodeID);
  		edgeList.push_back(e); 
    }

    void clearEdges(){
    	edgeList.clear();
    }

    float getRecoveryTime()
	{
		return recovery_time;
	}

	float UpdateRecoveryTime(float time)
	{
		recovery_time -= time;
	}
};

class Graph{

	vector<Node> nodes;

public:

	void addNode(Node newNode){
		bool check = checkIfNodeExistsByID(newNode.getID());
		if (check == true)
		{
			cout << "Node already exists" << endl;
			//break;
		}
		else
		{
			nodes.push_back(newNode);
		}
	}

	int NumberOfEdgesByNode(int currentNode){
		int temp2;

		temp2 = nodes.at(currentNode).edgeList.size();
		return temp2;
	}

	float DegreeOfNodeToAttach(int node_to_attach){
		return nodes.at(node_to_attach).edgeList.size();
	}

	//Removed weight from addEdgeID.
	void addEdgeByID(int fromNode, int toNode) {
		int CurrentEdgeNumber = 0;
    	bool check1 = checkIfNodeExistsByID(fromNode);
    	bool check2 = checkIfNodeExistsByID(toNode);

    	bool check3 = checkIfEdgeExistsByID(fromNode, toNode);
    	if ((check1 && check2 == true)) {

        	if (check3 == true) {
        		cout << "Edge between " << getNodeByID(fromNode).getID() << "(" << fromNode << ") and " << getNodeByID(toNode).getID() << "(" << toNode << ") Already Exist" << endl;
      		} 
      		else {
				for (int i = 0; i < nodes.size(); i++) {

            		if (nodes.at(i).getID() == fromNode) {
            			Edge e(toNode);
            			Edge e2(fromNode);
            			nodes.at(i).edgeList.push_back(e);
            			CurrentEdgeNumber = nodes.at(i).edgeList.size();
            			nodes.at(toNode).edgeList.push_back(e2);
          			} 
        		}
        	}
    	} 
    	else {
    		cout << "Invalid Node ID entered.";
    	}
    }

	bool checkIfNodeExistsByID(int nID){
		bool flag = false;

		for (int i = 0; i < nodes.size(); ++i)
		{
			if (nodes.at(i).getID()==nID)
			{
				return true;
			}
		}
		return flag;
	}

	bool checkIfEdgeExistsByID(int fromNode, int toNode) {
    	Node v = getNodeByID(fromNode);
    	list < Edge > e;
    	e = v.getEdgeList();
    	bool flag = false;
    	
    	for (auto it = e.begin(); it != e.end(); it++) {
      		if (it -> getDestinationID() == toNode) {
        		flag = true;
        		return flag;
        		break;
      		}

    	}

    	Node v2 = getNodeByID(toNode);
    	list < Edge > e2;
    	e2 = v2.getEdgeList();
    	
    	for (auto it = e2.begin(); it != e2.end(); it++) {
      		if (it -> getDestinationID() == fromNode) {
        		flag = true;
        		return flag;
        		break;
      		}

    	}
    	return flag;
  	}

	Node getNodeByID(int vid) {
		Node temp;
		for (int i = 0; i < nodes.size(); i++) {
			temp = nodes.at(i);
			if (temp.getID() == vid) {
			return temp;
			}
		}
		return temp;
	}

	void clearGraph(){
		for (int i = 0; i < nodes.size(); i++) {
		    Node temp;
		    temp = nodes.at(i);
		    temp.clearEdges();
		    
		}
		nodes.clear();
	}

	void printGraph() {
	    for (int i = 0; i < nodes.size(); i++) {
		    Node temp;
		    temp = nodes.at(i);
		    temp.printEdgeList();
		    temp.printNumberEdges();
	    }
    }

    void PrepareToRewire(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	temp.clearEdges();
    }

    void NodeInfection(int node)
    {
    	nodes.at(node).setStatus(1);
    }

    void NodeSusceptible(int node)
    {
    	nodes.at(node).setStatus(0);
    }

    void RecoveryTime(int node)
	{
		nodes.at(node).UpdateRecoveryTime(1);
	}
    
    vector<int> EdgesOfNode(int node){
        Node v = getNodeByID(node);
        list < Edge > e;
        e = v.getEdgeList();
        
        vector<int> arr;
        
        for (auto it = e.begin(); it != e.end(); it++) {
            arr.push_back(it -> getDestinationID());
        }
        
        return arr;
    }


    int checkStatus(int node)
    {
    	Node temp;
    	temp = nodes.at(node);
    	return temp.getStatus();
    }

    float checkRecoveryTime(int node)
	{
		return nodes.at(node).getRecoveryTime();
	}
    
};


int main(){

	Graph g;
	srand(time(0));

	int simulations = 1;
	int sim_number;
	int m = 10; // ki -> degree of old nodes... m -> degree of new node.
	int m0 = m+1, N = 1000;

	for (sim_number = 0; sim_number < simulations; sim_number++)
	{	
		g.clearGraph();

		int contador = 0;
		int sum_ki = 0;
		int new_nodes = m0;
		int random_node;
		double Pk = 0;
		float dice = 0;

		// Initial nodes
		for (int i = 0; i < m0; ++i)
		{
			Node v;
			v.setID(i);
			g.addNode(v);
		}
	
		// Now we need to connect them somehow. This is fully connected.
		for (int i = 0; i < m0; ++i)
		{
			for (int j = 0; j < m0; ++j)
			{
				if (i != j)
				{
					g.addEdgeByID(i,j); //(from, to)
				}
			}
		}
		//g.printGraph();

		// Add a new node
		while(contador != N)
		{
			// Add a new node
			Node v;
			v.setID(new_nodes);
			g.addNode(v);
			//---------------
			
			// Sum the edges of all existing nodes.
			sum_ki = 0;
			for (int i = 0; i < new_nodes; ++i)
			{
				sum_ki += g.NumberOfEdgesByNode(i); //sum_ki is twice the total edges.
			}

			int conections_per_new_node = 0;

			// Pick a random node, check the prob to connect and throw a random number, if its lower than the prob, connect it.
			while(conections_per_new_node < m)
			{
				random_node = rand() % new_nodes;
				Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1);
				
				dice = (float)rand()/RAND_MAX;
				
				if (Pk > dice && g.checkIfEdgeExistsByID(new_nodes,random_node)==false)
				{
					g.addEdgeByID(new_nodes,random_node);
					conections_per_new_node++;
				}
			}
			contador++;
			new_nodes++;
		}
		//g.printGraph();

		//----------------------------------------------------------------------------------------------------------------
		//                                         Estatistical Distribution (**)
		//----------------------------------------------------------------------------------------------------------------
		list<float> probabilidades_nodos;
		double total_prob = 0;
		double prob_nodo[N+5], prob_nodo_normal[N+5];
		double tabla[N+5];

		// Probability of choosing nodes based on the amount of edges.
		for (int i = 0; i < N; i++)
		{
			prob_nodo[i] = (double)1.0/g.NumberOfEdgesByNode(i);
			total_prob += prob_nodo[i];
		}

		// Normalized probability of nodes.
		for (int i = 0; i < N; i++)
		{
			prob_nodo_normal[i] = prob_nodo[i]/total_prob;
		}
		
		tabla[0] = prob_nodo_normal[0];

		// Generates the probability of each node inside the distribution.
		for (int i = 1; i < N; i++)
		{
			tabla[i] = tabla[i-1] + prob_nodo_normal[i];
		}
		
		//----------------------------------------------------------------------------------------------------------------
		//                                             Rewiring section
		//----------------------------------------------------------------------------------------------------------------
		int amount_of_nodes_to_rewire = 500;
		double choice;
		float rewiring_rate[N+1];
        int rewiring_node = 0, numberOfEdges = 0, contador2 = 0;
		
		// Generates a list containing the probability of choosing each node for rewiring.
		for (int i = 0; i < N; ++i)
		{
			rewiring_rate[i] = m/g.DegreeOfNodeToAttach(i);
		}
		
		// Rewire nodes until the desired amount.
		while(contador2 < amount_of_nodes_to_rewire)
		{
			bool flag = false;

			while(flag != true)
			{
				rewiring_node = rand() % N;
				choice = (float)rand()/RAND_MAX;
				// Removes the edges from the current selected node.
				if(choice < rewiring_rate[rewiring_node])
				{
					numberOfEdges = g.NumberOfEdgesByNode(rewiring_node);
					g.PrepareToRewire(rewiring_node);
					flag = true;
					contador2++;
				}
			}

			int conections_per_rewired_node = 0;

			// Pick a random node, check the prob to connect and throw a random number, if its lower than the prob, connect it.
			while(conections_per_rewired_node < numberOfEdges)
			{
				random_node = (rand() % (N-1));
				if (rewiring_node != random_node)
				{
					Pk = pow((g.DegreeOfNodeToAttach(random_node)/sum_ki),1); 

					dice = (float)rand()/RAND_MAX;
			
					if (Pk > dice && g.checkIfEdgeExistsByID(rewiring_node,random_node)==false)
					{
						g.addEdgeByID(rewiring_node,random_node);
						conections_per_rewired_node++;
					}
				}
			}
		}

		//----------------------------------------------------------------------------------------------------------------
		//                  This section is designed to plot the characteristic curve of the SF network.
		//----------------------------------------------------------------------------------------------------------------
		/*
		int ListNumberOfEdges[N+1];
		int FreqDist[N+1];

		// We need to initialize the array to avoid getting random elements from accupied memory spaces.
		for (int i = 0; i < N+2; ++i)
		{
			FreqDist[i] = 0;
		}

		for (int i = 0; i < N+1; ++i)
		{
			ListNumberOfEdges[i] = g.NumberOfEdgesByNode(i); // Node "i" has "g.NumberOfEdgesByNode(i)" Edges.
			FreqDist[ListNumberOfEdges[i]]++; //Now ListNumberOfEdges[i] is the index of the FreqDist array and we add one to that index every time we find it.
		}

		//----------------------------------------------------------------------------------------------------------------
		//                                             Printing section
		//----------------------------------------------------------------------------------------------------------------
		string output = "sf2_network_" + to_string(N) + "_rewired_" + to_string(amount_of_nodes_to_rewire) + "_" + to_string(sim_number+1);
		ofstream sf(output);

		for (int i = 0; i < N+1; ++i)
		{
			if (i >= m)
			{
				sf << i << " " << FreqDist[i] << endl;
			}
		}
		*/
		//----------------------------------------------------------------------------------------------------------------
		
		//sf.close();
	
	float odds, beta = 1.0;
	list<int> infecteds;
	int initial_infecteds = 0;
	
	// We select the amount of initial infections in the network
	do
	{
		int random_node = (rand() % (N-1));
		float rn_probability = rewiring_rate[random_node];
		odds = (float)rand()/RAND_MAX;

		// Infect the node if it's not already infected. Here beta is not relevant because this is the initial infection.
		if (rn_probability > odds && g.checkStatus(random_node) != 1)
		{
			g.NodeInfection(random_node);
			initial_infecteds++;
			cout << "node number " << random_node << " has been initially infected" << endl;
		}
	} while (initial_infecteds < 1);
	

	int betastr = beta * 10; // This was just used to make my life easier while naming files.
    string output2 = "SIS_sf_network_" + to_string(N) + "_rewired_" + to_string(amount_of_nodes_to_rewire) + "_m10_sim_" + to_string(sim_number) + "_beta_" + to_string(betastr);
    ofstream sfsir(output2);
	
	int status_for_switch;
	int t = 0, tf = 300;
	int susceptibles = 0, infkts = 0;

			// Count the amount of susceptible and infected nodes in the network.
			for (int i = 0; i < N; i++)
			{
				if (g.checkStatus(i) == 0)
				{
					susceptibles++;
				}
				else if (g.checkStatus(i) == 1)
				{
					infkts++;
				}
			}
		sfsir << t << " " << susceptibles << " " << infkts << endl; 

	// Disease spread section. Every Monte Carlo step is equivalent to N actions.
	// This is, checking the contact between a node and its edges for N nodes
	// WHich are chosen using the statistical probability defined in (**).

    while(t < tf){
		// We select a node from 
		for (int i = 0; i < N; i++)
		{
			double num_aleatorio = (double)rand()/RAND_MAX;
			int selected_node = 0;

			for (int j = 1; j <= N; j++)
			{
				if (tabla[j-1] < num_aleatorio <= tabla[j])
				{
					selected_node = j;
					break;
				}
				
			}
			vector<int> Connections = g.EdgesOfNode(selected_node); // Saving the edges from this node.
			status_for_switch = g.checkStatus(selected_node);

			switch(status_for_switch)
			{
				case 0:
				for (int j = 0; j < Connections.size(); ++j)
				{
					if (g.checkStatus(Connections[j]) == 1)
					{
						odds = (float)rand()/RAND_MAX;
							
						if(odds < beta){
							g.NodeInfection(selected_node); //supposing the chances of getting the infection are X% when in contact.
						}
					}
				}
				break;
				case 1:
				g.RecoveryTime(selected_node);

				if(g.checkRecoveryTime(selected_node) <= 0){
					g.NodeSusceptible(selected_node);
				}
				break;
			}
			Connections.clear();
		}

        //====================================================  PRINT SIS  ===================================================
		susceptibles = 0; 
		infkts = 0;

		for (int i = 0; i < N; i++)
		{
			if (g.checkStatus(i) == 0)
			{
				susceptibles++;
			}
			else if (g.checkStatus(i) == 1)
			{
				infkts++;
			}
		}
		
        sfsir << t << " " << susceptibles << " " << infkts << endl;
        
        
        //====================================================================================================================
        t++;
    }
	}

	return 0;
}//end
