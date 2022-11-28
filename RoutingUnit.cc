/*
 * Copyright (c) 2008 Princeton University
 * Copyright (c) 2016 Georgia Institute of Technology
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Authors: Niket Agarwal
 *          Tushar Krishna
 */


// #include "mem/ruby/network/garnet2.0/RoutingUnit.hh"

// #include "base/cast.hh"
// #include "base/logging.hh"
// #include "mem/ruby/network/garnet2.0/InputUnit.hh"
// #include "mem/ruby/network/garnet2.0/Router.hh"
// #include "mem/ruby/slicc_interface/Message.hh"


#include "mem/ruby/network/garnet2.0/RoutingUnit.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "base/cast.hh"
#include "debug/RubyNetwork.hh"
#include "mem/ruby/network/garnet2.0/InputUnit.hh"
#include "mem/ruby/network/garnet2.0/OutputUnit.hh"
#include "mem/ruby/network/garnet2.0/Router.hh"
#include "mem/ruby/slicc_interface/Message.hh"
#include "mem/ruby/network/garnet2.0/flit.hh"
#include "cpu/kvm/base.hh"
#include <Python.h>
#include "pyhelper.hpp"
#include <unistd.h>
#include <math.h>



RoutingUnit::RoutingUnit(Router *router)
{
    m_router = router;
    m_routing_table.clear();
    m_weight_table.clear();
}

void
RoutingUnit::addRoute(const NetDest& routing_table_entry)
{
    m_routing_table.push_back(routing_table_entry);
}

void
RoutingUnit::addWeight(int link_weight)
{
    m_weight_table.push_back(link_weight);
}

/*
 * This is the default routing algorithm in garnet.
 * The routing table is populated during topology creation.
 * Routes can be biased via weight assignments in the topology file.
 * Correct weight assignments are critical to provide deadlock avoidance.
 */

int
RoutingUnit::lookupRoutingTable(int vnet, NetDest msg_destination)
{
    // First find all possible output link candidates
    // For ordered vnet, just choose the first
    // (to make sure different packets don't choose different routes)
    // For unordered vnet, randomly choose any of the links
    // To have a strict ordering between links, they should be given
    // different weights in the topology file

    int output_link = -1;
    int min_weight = INFINITE_;
    std::vector<int> output_link_candidates;
    int num_candidates = 0;

    // Identify the minimum weight among the candidate output links
    for (int link = 0; link < m_routing_table.size(); link++) {
        if (msg_destination.intersectionIsNotEmpty(m_routing_table[link])) {

        if (m_weight_table[link] <= min_weight)
            min_weight = m_weight_table[link];
        }
    }

    // Collect all candidate output links with this minimum weight
    for (int link = 0; link < m_routing_table.size(); link++) {
        if (msg_destination.intersectionIsNotEmpty(m_routing_table[link])) {

            if (m_weight_table[link] == min_weight) {

                num_candidates++;
                output_link_candidates.push_back(link);
            }
        }
    }

    if (output_link_candidates.size() == 0) {
        fatal("Fatal Error:: No Route exists from this Router.");
        exit(0);
    }

    // Randomly select any candidate output link
    int candidate = 0;
    if (!(m_router->get_net_ptr())->isVNetOrdered(vnet))
        candidate = rand() % num_candidates;

    output_link = output_link_candidates.at(candidate);
    return output_link;
}


void
RoutingUnit::addInDirection(PortDirection inport_dirn, int inport_idx)
{
    m_inports_dirn2idx[inport_dirn] = inport_idx;
    m_inports_idx2dirn[inport_idx]  = inport_dirn;
}

void
RoutingUnit::addOutDirection(PortDirection outport_dirn, int outport_idx)
{
    m_outports_dirn2idx[outport_dirn] = outport_idx;
    m_outports_idx2dirn[outport_idx]  = outport_dirn;
}

// outportCompute() is called by the InputUnit
// It calls the routing table by default.
// A template for adaptive topology-specific routing algorithm
// implementations using port directions rather than a static routing
// table is provided here.

int
RoutingUnit::outportCompute(flit *t_flit, RouteInfo route, int inport,
                            PortDirection inport_dirn)
{
    int outport = -1;

    if (route.dest_router == m_router->get_id()) {

        // Multiple NIs may be connected to this router,
        // all with output port direction = "Local"
        // Get exact outport id from table
        outport = lookupRoutingTable(route.vnet, route.net_dest);
        return outport;
    }

    // Routing Algorithm set in GarnetNetwork.py
    // Can be over-ridden from command line using --routing-algorithm = 1
    RoutingAlgorithm routing_algorithm =
        (RoutingAlgorithm) m_router->get_net_ptr()->getRoutingAlgorithm();
	std::set<int> paths;
    switch (routing_algorithm) {
        case TABLE_:  outport =
            lookupRoutingTable(route.vnet, route.net_dest); break;
        case XY_:     outport =
            outportComputeXY(route, inport, inport_dirn); break;
        // any custom algorithm
        case OE_: outport =
            outportComputeOE(route, inport, inport_dirn, paths); break;
    case Q_R_: outport = outportComputeQ_Routing(t_flit, inport, inport_dirn); break;
        default: outport =
            lookupRoutingTable(route.vnet, route.net_dest); break;
    }

    assert(outport != -1);
    return outport;
}

// XY routing implemented using port directions
// Only for reference purpose in a Mesh
// By default Garnet uses the routing table
int
RoutingUnit::outportComputeXY(RouteInfo route,
                              int inport,
                              PortDirection inport_dirn)
{
    PortDirection outport_dirn = "Unknown";

    int M5_VAR_USED num_rows = m_router->get_net_ptr()->getNumRows();
    int num_cols = m_router->get_net_ptr()->getNumCols();
    assert(num_rows > 0 && num_cols > 0);

    int my_id = m_router->get_id();
    int my_x = my_id % num_cols;
    int my_y = my_id / num_cols;

    int dest_id = route.dest_router;
    int dest_x = dest_id % num_cols;
    int dest_y = dest_id / num_cols;

    int x_hops = abs(dest_x - my_x);
    int y_hops = abs(dest_y - my_y);

    bool x_dirn = (dest_x >= my_x);
    bool y_dirn = (dest_y >= my_y);

    // already checked that in outportCompute() function
    assert(!(x_hops == 0 && y_hops == 0));

    if (x_hops > 0) {
        if (x_dirn) {
            assert(inport_dirn == "Local" || inport_dirn == "West");
            outport_dirn = "East";
        } else {
            assert(inport_dirn == "Local" || inport_dirn == "East");
            outport_dirn = "West";
        }
    } else if (y_hops > 0) {
        if (y_dirn) {
            // "Local" or "South" or "West" or "East"
            assert(inport_dirn != "North");
            outport_dirn = "North";
        } else {
            // "Local" or "North" or "West" or "East"
            assert(inport_dirn != "South");
            outport_dirn = "South";
        }
    } else {
        // x_hops == 0 and y_hops == 0
        // this is not possible
        // already checked that in outportCompute() function
        panic("x_hops == y_hops == 0");
    }

    return m_outports_dirn2idx[outport_dirn];
}

// Template for implementing custom routing algorithm
// using port directions. (Example adaptive)
// int
// RoutingUnit::outportComputeCustom(RouteInfo route,
//                                  int inport,
//                                  PortDirection inport_dirn)
// {
//   PortDirection outport_dirn = "Unknown";
    
//     int num_rows = m_router->get_net_ptr()->getNumRows();
//     int num_cols = m_router->get_net_ptr()->getNumCols();
//     int num_tiles = num_rows * num_cols;
//     assert(num_rows > 0 && num_cols > 0);

//     int id = m_router->get_id();
//     int source_id = route.src_router; //check this 
//     int dest_id = route.dest_router;  
//     // id = 11 , cur_x = 2 , cur_y = 3
//     int cur_yco = id / num_cols;
//     int cur_xco = id % num_cols;
    
//     //int src_yco = source_id / num_cols;
//     int src_xco = source_id / num_rows;
    
//     int dst_yco = dest_id / num_cols;
//     int dst_xco = dest_id % num_cols;

//     int dif_yco = dst_yco - cur_yco;
//     int dif_xco = dst_xco - cur_xco;

//     if (dif_yco == 0 && dif_xco == 0)   // cur node is dest
//         return m_outports_dirn2idx[outport_dirn];

//     if (dif_xco == 0)   // cur node in same col as dest
//     {
//         if (dif_yco < 0)
//         {
//             if (inport_dirn != "South" && !borderS(id, num_tiles, num_cols))        // 180-degree turn prohibited
//                 outport_dirn = "South"; // allow to route s
//         }
//         else
//         {
//             if (inport_dirn != "North" && !borderN(id, num_cols))
//                 outport_dirn = "North"; // allow to route N
//         }
//     }
//     else
//     {
//         if (dif_xco > 0)    // E-bound pkt
//         {
//             if (dif_yco == 0)   // cur in same row as dest
//             {
//                 if (inport_dirn != "East")
//                     outport_dirn = "East";
//             }
//             else
//             {
//                 if (cur_xco % 2 != 0 || cur_xco == src_xco) // N/S turn allowed only in odd col.
//                 {
//                     if (dif_yco < 0 && !borderS(id, num_tiles, num_cols) && inport_dirn != "South")
//                         outport_dirn = "South";
//                     else if (!borderN(id, num_cols) && inport_dirn != "North")
//                         outport_dirn = "North";
//                 }                
//                 if (dst_xco % 2 != 0 || dif_xco != 1)   // allow to go E only if dest is odd col
//                 {
//                     if (inport_dirn != "East")
//                         outport_dirn = "East";  // because N/S turn not allowed in even col.
//                 }
//             }
//         }
//         else    // W-bound
//         {
// 	  if ((cur_xco % 2) == 0) {  // allow to go N/S only in even col. because N->W and S->W not allowed in odd col.
//                 if (dif_yco <= 0 && !borderS(id, num_tiles, num_cols) && inport_dirn != "South") // = 0 to allow non minimal path
//                     outport_dirn = "South";
//                 else if (!borderN(id, num_cols) && inport_dirn != "North")
//                     outport_dirn = "North";
// 	  }
//         }
//     }
//     return m_outports_dirn2idx[outport_dirn];
//   //panic("%s placeholder executed", __FUNCTION__);
// }

int RoutingUnit::outportComputeOE(RouteInfo route,
                              int inport,
                              PortDirection inport_dirn, std::set<int> &paths)
{
    ////std::cout<<"In outportComputeOE"<<std::endl;
    PortDirection outport_dirn = "Unknown";
    
    int M5_VAR_USED num_rows = m_router->get_net_ptr()->getNumRows();
    int num_cols = m_router->get_net_ptr()->getNumCols();
    assert(num_rows > 0 && num_cols > 0);

    int my_id = m_router->get_id();
    int my_x = my_id % num_cols;
    int my_y = my_id / num_cols;

    int dest_id = route.dest_router;
    int dest_x = dest_id % num_cols;
    int dest_y = dest_id / num_cols;

    int src_id = route.src_router;
    int src_x = src_id % num_cols;
    //    int src_y = src_id / num_cols;

    int x_hops = dest_x - my_x;
    int y_hops = dest_y - my_y;

    // already checked that in outportCompute() function
    assert(!(x_hops == 0 && y_hops == 0));
	
	if(x_hops == 0) {
		if(y_hops < 0) {
			outport_dirn = "South";
			paths.insert(2);
		}
		else {
			outport_dirn = "North";
			paths.insert(0);
		}
	}
	else {
		if(x_hops > 0) {
			if(y_hops == 0) {
				outport_dirn = "East";
				paths.insert(1);
			}
			else {
				if(my_x % 2 != 0 || my_x == src_x) {
					if(y_hops < 0) {
						outport_dirn = "South";
						paths.insert(2);
					}
					else {
						outport_dirn = "North";
						paths.insert(0);
					}
				}
				if(dest_x % 2 != 0 || x_hops != 1) {
					outport_dirn = "East";
					paths.insert(1);
				}
			}
		}
		else {
			outport_dirn = "West";
			paths.insert(3);
			if(y_hops != 0) {
				if(my_x % 2 == 0) {
					if(y_hops < 0) {
						outport_dirn = "South";
						paths.insert(2);
					}
					else {
						outport_dirn = "North";
						paths.insert(0);
					}
				}
			}	
		}
	}
    return m_outports_dirn2idx[outport_dirn];
}


bool RoutingUnit::borderN(int id, int num_cols) {
    return id<num_cols;
}
bool RoutingUnit::borderS(int id, int num_tiles, int num_cols) {
    return id>= num_tiles - num_cols;
}

bool RoutingUnit::borderE(int id, int num_cols) {
    return ((id+1)% num_cols ==0);
}

bool RoutingUnit::borderW(int id, int num_cols) {
    return ((id%num_cols)==0);
}

bool RoutingUnit::bored(int id, int num_cols, int num_tiles) {
    return borderN(id, num_cols) || borderS(id, num_tiles, num_cols) || borderE(id, num_cols) || borderW(id,num_cols) ;
}

#define EPSILON 0.3
#define epsilon 0.3
#define NROUTERS 64
#define NACTIONS 4
#define GRIDSIZE 8
#define LEARNINGRATE 0.5
#define DISCOUNTRATE 0.7
#define EPSILON_DECAY 0.99975
#define MIN_EPSILON 0.001

int RoutingUnit::epsilon_greedy(std::vector<std::vector<std::vector<double>>> Q, int state, int destination) {
        ////std::cout<<"In epsilon_greedy size of q : "<<Q[state][destination].size()<<" q values : "<<Q[state][destination][0]<<" "<<Q[state][destination][1]<<" "<<Q[state][destination][2]<<" "<<Q[state][destination][3]<<" "<<std::endl;
	float p = (float) rand() / RAND_MAX;
	if(p > EPSILON) {
		int optimalAction = std::distance(Q[state][destination].begin(), std::min_element(Q[state][destination].begin(), Q[state][destination].end()));
		//std::cout<<"state: "<<state<<" destination : "<<destination<<" Optimal Action: "<<optimalAction<<std::endl;
		return optimalAction;
	}

        int randomAction = rand() % 4;
//	randomAction = temp % 4;
	////std::cout<<"Rchd RandomAction: "<<randomAction<<std::endl;
	return randomAction;
}

int RoutingUnit::get_action(std::vector<std::vector<std::vector<double>>> Q, int state, int destination, std::set<int> paths) {
	//std::cout<<"In get_action size of paths is "<<paths.size();
	double min_q = std::numeric_limits<double>::max();
	int min_state = -1;
	for(auto & path:paths) {
		//std::cout<<" path: "<<path<<"  Q val: "<<Q[state][destination][path];
		if(Q[state][destination][path] <=min_q) {
			min_q = Q[state][destination][path];
			min_state = path;
		}
	}
	//std::cout<<" min_state: "<<min_state<<" min_q: "<<min_q<<std::endl;

	return min_state;
}


int RoutingUnit::outportComputeQ_Routing(flit *t_flit, int inport, PortDirection inport_dirn) {
	Tick src_queueing_delay = t_flit->get_src_delay();
    Tick dest_queueing_delay = (curTick() - t_flit->get_dequeue_time());
    Tick queueing_delay = src_queueing_delay + dest_queueing_delay;
//	//std::cout<<"CURRENT TICK "<<curTick()<<"\n";	
	static int SEED = 0;
	if(SEED == 0){
		srand(time(NULL));
//		//std::cout << "Reaching srand()" << std::endl;
		SEED += 1;
	}
	static int iter = 0;
	////std::cout << "Number of iterations: " << iter << std::endl;
	iter++;
	
    PortDirection outport_dirn = "Unknown";
    static bool isQTableInitialized = false;
	static std::vector<std::vector<std::vector<double>>> Q(NROUTERS, std::vector<std::vector<double>>(NROUTERS, std::vector<double> (NACTIONS, INT_MAX)));
	if(!isQTableInitialized){
		isQTableInitialized = true;
		if(access("Q_Table.txt", F_OK) == 0) {
			std::ifstream f_qTable {"Q_Table.txt"};
			for(int i=0;i<NROUTERS;++i) {
				for(int j=0;j<NROUTERS;++j) {
					for(int k=0;k<NACTIONS;++k) {
						f_qTable >> Q[i][j][k];
					}
				}
			}		
		}
	}

	//---Initializing Q-Table---
	RouteInfo route = t_flit->get_route();

	//---Geting source and destination router details
	int M5_VAR_USED num_rows = m_router->get_net_ptr()->getNumRows();
	int num_cols = m_router->get_net_ptr()->getNumCols();
    assert(num_rows > 0 && num_cols > 0);
    int my_id = m_router->get_id();
    int my_x = my_id % num_cols;
    int my_y = my_id / num_cols;

    int dest_id = route.dest_router;

    int src_id = route.src_router;
	
	//int action = epsilon_greedy(Q, my_id, dest_id);//finding the minimum q value between source and destination 
	int prev_router_id;

	int temp_x = 0;
	int temp_y = 0;
	if(inport_dirn == "North"){
		if(my_y < num_rows - 1) {
			temp_y = my_y + 1;
			temp_x = my_x;
		}
	}
	else if(inport_dirn == "South") {
		if(my_y > 0) {
			temp_y = my_y - 1;
			temp_x = my_x;
		}
	
	}
	else if(inport_dirn == "East") {
		if(my_x < num_cols -1) {
			temp_x = my_x + 1;
			temp_y = my_y;
		}
	}
	else if(inport_dirn == "West"){
		if(my_x > 0) {
			temp_x = my_x - 1;
			temp_y = my_y;
		}
	
	}
	else{
	}

	prev_router_id = temp_y * num_cols + temp_x;

	//int my_action = get_output_port(Q, my_id, dest_id);

	std::set<int> paths;

	outportComputeOE(route, inport, inport_dirn, paths);

	int my_action = get_action(Q, my_id, dest_id, paths);

	//std::cout<<" action : "<<action <<" my_action: "<<my_action<<std::endl;

	if(my_action == 0 && my_y < num_rows-1) {
		outport_dirn = "North";
	}
	else if(my_action == 1 && my_x < num_cols-1) {
		outport_dirn = "East";
	}
	else if(my_action == 2 && my_y>0) {
		outport_dirn = "South";
	}
	else if(my_action == 3 && my_x>0){
		outport_dirn = "West";
	}
	else{
		long long random = rand();
		my_action = random % 4;
	}		
	
	if(my_id == src_id) {
		return m_outports_dirn2idx[outport_dirn];
	}


	int next_id = m_outports_dirn2idx[outport_dirn];

	int link_latency = find_link_latency(my_id, next_id);
	
	int prev_action = 0;
	if(inport_dirn == "North") {
		prev_action = 2;
	}
	else if(inport_dirn == "East") {
		prev_action = 3;
	}
	else if(inport_dirn == "South") {
		prev_action = 0;
	}
	else {
		prev_action = 1;
	}
	
	// //std::cout << "Iterations: " << iter <<std::endl;
	//Updating Q-Table
	////std::cout << "Queueing Delay: " << queueing_delay << std::endl;
	double Qy_min = *std::min_element(Q[my_id][dest_id].begin(), Q[my_id][dest_id].end());
	//std::cout << "Q-tavle value before updation: " << Q[prev_router_id][dest_id][prev_action] << std::endl;
	// //std::cout << "Qy_min: " << Qy_min << " Queueing delay: " << queueing_delay << std::endl;
	Q[prev_router_id][dest_id][prev_action] = Q[prev_router_id][dest_id][prev_action] + LEARNINGRATE * (DISCOUNTRATE*Qy_min + ((int)(queueing_delay) + link_latency) - Q[prev_router_id][dest_id][prev_action]);
	//Q[prev_router_id][dest_id][prev_action] = Q[prev_router_id][dest_id][prev_action] + LEARNINGRATE * (Qy_min + ((int)(queueing_delay) + link_latency) - Q[prev_router_id][dest_id][prev_action]);
	if(my_id == dest_id){
	  Q[prev_router_id][dest_id][prev_action] = Q[prev_router_id][dest_id][prev_action] + LEARNINGRATE * (((int)(queueing_delay)+link_latency) - Q[prev_router_id][dest_id][prev_action]);
		int outport = lookupRoutingTable(route.vnet, route.net_dest);
	    return outport;
	}
	static bool file_dumped = false;
	//std::cout<<"curTick "<<curTick()<<std::endl;
	if(curTick() >2000) {			
		if(!file_dumped) {
			file_dumped = true;	
			std::cout<<"Current tick is "<< curTick()<<std::endl;
			std::ofstream f_qTable;
			////std::cout<<"Updating File\n";
			f_qTable.open("Q_Table.txt", std::ios::out | std::ios::trunc);
			for(int i=0;i<NROUTERS;++i) {
				for(int j=0;j<NROUTERS;++j) {
					for(int k=0;k<NACTIONS;++k) {
						f_qTable << Q[i][j][k] << " ";
					}
					f_qTable << "\n";
				}
				f_qTable << "\n";
			}
		}	
	}
    return m_outports_dirn2idx[outport_dirn];
}


int RoutingUnit::find_link_latency(int src, int dest) {
  if((src_nodes_1.find(src)!= src_nodes_1.end()) && (dest_nodes_1.find(dest)!= dest_nodes_1.end()))
    return 500;
  else if ((src_nodes_2.find(src)!= src_nodes_2.end()) && (dest_nodes_2.find(dest)!= dest_nodes_2.end()))
    return 1000;
  else if ((src_nodes_3.find(src)!= src_nodes_3.end()) && (dest_nodes_3.find(dest)!= dest_nodes_3.end()))
    return 1500;
  else
    return 1;
}

