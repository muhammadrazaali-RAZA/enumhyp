#include "hypergraph.h"

#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <mpi.h>
#include <thread>
#include <vector>
#include <queue>    // For std::queue
#include <tuple>    // For std::tuple
#include <unordered_map>

#define NOT_EXTENDABLE 0
#define EXTENDABLE 1
#define MINIMAL 2

Hypergraph::Hypergraph() {
}

Hypergraph::Hypergraph(const Hypergraph &other) {
	m_edges = other.m_edges;
	m_num_vertices = other.m_num_vertices;
}

Hypergraph::Hypergraph(std::string path) {
	m_edges.clear();
	std::ifstream infile(path);
	std::string s_num_vertices;
	getline(infile, s_num_vertices);
	int num_vertices = std::stoi(s_num_vertices);
	if (num_vertices <= 0) throw "Tried to read graph with less than one vertex from file!";
	m_num_vertices = num_vertices;
	edge_set edges;
	for (std::string line; getline(infile, line); ) {
		std::stringstream ss(line);
		std::string node;
		edge e((edge::size_type) num_vertices);
		while (std::getline(ss, node, ',')) {
			int i_vertex = std::stoi(node);
			if (i_vertex < 0 || i_vertex >= num_vertices) throw "Found corrupt edge while reading graph from file!";
			e[(edge::size_type) i_vertex] = 1;
		}
		edges.insert(e);
	}
	for (auto e : edges) m_edges.push_back(e);
}

Hypergraph::Hypergraph(const Table &t) {
	if (t.empty()) return;
	m_num_vertices = (int)t.m_records[0].size();
	m_edges = t.edges();
	minimize();
}

Hypergraph::Hypergraph(int num_vertices, edge_vec edges) {
	m_num_vertices = num_vertices;
	m_edges = edges;
}

Hypergraph::~Hypergraph() {
}

void Hypergraph::print_edges() const {
	print_edge_vec(m_edges);
}

void Hypergraph::save(std::string path) const {
	std::ofstream outfile;
	outfile.open(path);
	outfile << m_num_vertices << std::endl;
	for (auto e : m_edges) {
		auto num_vertices = e.count();
		if (num_vertices == 0) {
			std::cerr << "Graph contains empty edge!" << std::endl;
			continue;
		}
		auto index = e.find_first();
		edge::size_type i = 0;
		while (i < num_vertices - 1) {
			outfile << index << ",";
			index = e.find_next(index);
			++i;
		}
		outfile << index << std::endl;
	}
	outfile.close();
}

bool Hypergraph::is_hitting_set(const edge &h) const {
	for (auto e : m_edges) if (!e.intersects(h)) return false;
	return true;
}

int Hypergraph::extendable(const edge &x, const edge &y) {
	if (m_configuration.collect_oracle_statistics) m_oracle_timestamp = Clock::now();
	// 2
	if (x.none()) {
		// 3
		if (is_hitting_set(~y)) {
			if (m_configuration.collect_oracle_statistics) {
				auto now = Clock::now();
				m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "3", ns_string(m_oracle_timestamp, now), "", "", "", "", "", "", "", "" });
			}
			return EXTENDABLE;
		}
		// 4
		if (m_configuration.collect_oracle_statistics) {
			auto now = Clock::now();
			m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "4", ns_string(m_oracle_timestamp, now), "", "", "", "", "", "", "", "" });
		}
		return NOT_EXTENDABLE;
	}
	// 5
	edge_vec t;
	// 6
	std::vector<edge_vec> s(x.count(), edge_vec());
	std::vector<std::vector<edge_vec>::size_type> x_index_to_s_index(m_num_vertices, -1);
	std::vector<edge_vec>::size_type s_index = -1;
	for (auto x_index = x.find_first(); x_index != edge::npos; x_index = x.find_next(x_index)) x_index_to_s_index[x_index] = ++s_index;
	// 7
	for (edge e : m_edges) {
		edge intersection = e & x;
		// 9
		if (intersection.none()) {
			t.push_back(e - y);
			continue;
		}
		// 8
		if (intersection.count() == 1) s[x_index_to_s_index[intersection.find_first()]].push_back(e - y);
	}
	// 10
	for (auto sx : s) if (sx.empty()) {
		if (m_configuration.collect_oracle_statistics) {
			auto now = Clock::now();
			m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "10", ns_string(m_oracle_timestamp, now), "", "", "", "0", "0", "0", std::to_string(t.size()), std::to_string(total_number_of_vertices_in_t(t)) });
		}
		return NOT_EXTENDABLE;
	}
	// 11
	if (t.empty()) {
		if (m_configuration.collect_oracle_statistics) {
			auto now = Clock::now();
			m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "11", ns_string(m_oracle_timestamp, now), "", "", std::to_string(maximum_iteration_count(s)), std::to_string(s.size()), std::to_string(summed_sx_sizes(s)), std::to_string(total_number_of_vertices_in_s(s)), "0", "0" });
		}
		return MINIMAL;
	}
	// 12
	if (m_configuration.collect_oracle_statistics) {
		m_iteration_count = 0;
		m_oracle_bf_timestamp = Clock::now();
	}
	std::vector<edge::size_type> iteration_position(x.count(), 0);
	while (true) {
		if (m_configuration.collect_oracle_statistics) m_iteration_count++;
		// 13
		edge w(m_num_vertices);
		bool increase_next = true;
		for (std::vector<edge_vec>::size_type i_s = 0; i_s < s.size(); ++i_s) {
			w |= s[i_s][iteration_position[i_s]];
			if (increase_next) {
				++iteration_position[i_s];
				if (iteration_position[i_s] == s[i_s].size()) iteration_position[i_s] = 0;
				else increase_next = false;
			}
		}
		// 14
		bool all_no_subset = true;
		for (edge e : t) {
			if (e.is_subset_of(w)) {
				all_no_subset = false;
				break;
			}
		}
		if (all_no_subset) {
			if (m_configuration.collect_oracle_statistics) {
				auto now = Clock::now();
				m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "14", ns_string(m_oracle_timestamp, now), ns_string(m_oracle_bf_timestamp, now), std::to_string(m_iteration_count), std::to_string(maximum_iteration_count(s)), std::to_string(s.size()), std::to_string(summed_sx_sizes(s)), std::to_string(total_number_of_vertices_in_s(s)), std::to_string(t.size()), std::to_string(total_number_of_vertices_in_t(t)) });
			}
			return EXTENDABLE;
		}
		if (increase_next) break;
	}
	// 15
	if (m_configuration.collect_oracle_statistics) {
		auto now = Clock::now();
		m_oracle_stats.add_record({ edge_to_string(x), edge_to_string(y), "15", ns_string(m_oracle_timestamp, now), ns_string(m_oracle_bf_timestamp, now), std::to_string(m_iteration_count), std::to_string(maximum_iteration_count(s)), std::to_string(s.size()), std::to_string(summed_sx_sizes(s)), std::to_string(total_number_of_vertices_in_s(s)), std::to_string(t.size()), std::to_string(total_number_of_vertices_in_t(t)) });
	}
	return NOT_EXTENDABLE;
}

Hypergraph Hypergraph::enumerate(enumerate_configuration configuration) {
	m_configuration = configuration;
	if (m_configuration.collect_hitting_set_statistics) {
		m_hitting_set_stats.clear();
		m_hitting_set_stats.add_record({ "minimal_hitting_set", "delay_ns" });
	}
	if (m_configuration.collect_oracle_statistics) {
		m_oracle_stats.clear();
		m_oracle_stats.add_record({ "x", "y", "return_line", "total_time_ns", "bf_time_ns", "actual_iteration_count", "maximum_iteration_count", "s_size", "summed_sx_sizes", "total_number_of_vertices_in_s", "t_size", "total_number_of_vertices_in_t" });
	}
	Hypergraph h;
	if (m_configuration.implementation == "standard") h = Hypergraph(m_num_vertices, enumerate());
	else if (m_configuration.implementation == "legacy") h = Hypergraph(m_num_vertices, enumerate_legacy());
	else if (m_configuration.implementation == "brute_force") h = Hypergraph(m_num_vertices, brute_force_mhs());
	else std::cerr << "Implementation " << m_configuration.implementation << " not found!";
	save_statistics();
	return h;
}

edge_vec Hypergraph::enumerate() {
	edge_vec minimal_hitting_sets;
	if (m_configuration.collect_hitting_set_statistics) m_hitting_set_timestamp = Clock::now();
	enumerate(edge(m_num_vertices), edge(m_num_vertices), 0, minimal_hitting_sets);
	return minimal_hitting_sets;
}

void Hypergraph::enumerate(const edge &x, const edge &y, edge::size_type r, edge_vec &minimal_hitting_sets) {
	edge xv = x;
	xv[r] = 1;
	switch (extendable(xv, y)) {
	case MINIMAL:
	{
		minimal_hitting_sets.push_back(xv);
		if (m_configuration.collect_hitting_set_statistics) {
			auto now = Clock::now();
			m_hitting_set_stats.add_record({ edge_to_string(xv), ns_string(m_hitting_set_timestamp, now) });
			m_hitting_set_timestamp = Clock::now();
		}
		break;
	}
	case EXTENDABLE:
		enumerate(xv, y, r + 1, minimal_hitting_sets);
		break;
	case NOT_EXTENDABLE:
		edge yv = y;
		yv[r] = 1;
		enumerate(x, yv, r + 1, minimal_hitting_sets);
		return;
	}
	edge yv = y;
	yv[r] = 1;
	if (extendable(x, yv)) enumerate(x, yv, r + 1, minimal_hitting_sets);
}

edge_vec Hypergraph::enumerate_legacy() {
	edge_vec minimal_hitting_sets;
	if (m_configuration.collect_hitting_set_statistics) m_hitting_set_timestamp = Clock::now();
	enumerate_legacy(edge(m_num_vertices), edge(m_num_vertices), 0, minimal_hitting_sets);
	return minimal_hitting_sets;
}

void Hypergraph::enumerate_legacy(const edge &x, const edge &y, edge::size_type r, edge_vec &minimal_hitting_sets) {
    if (r == m_num_vertices) {
        minimal_hitting_sets.push_back(x);
        if (m_configuration.collect_hitting_set_statistics) {
            auto now = Clock::now();
            m_hitting_set_stats.add_record({ edge_to_string(x), ns_string(m_hitting_set_timestamp, now) });
            m_hitting_set_timestamp = Clock::now();
        }
        return;
    }

    edge xv = x;
    edge yv = y;
    xv[r] = 1;
    yv[r] = 1;

    assign_extendable_task(xv, y, x, yv,r+1,minimal_hitting_sets);
	return;
}


struct WorkData {
    std::vector<int> x_vec;
    std::vector<int> y_vec;
    int hitting_sets_size;
    std::vector<int> serialized_hitting_sets;
    unsigned long r;
};


void Hypergraph::assign_extendable_task(
    const edge &xv, const edge &y, const edge &x, const edge &yv, edge::size_type r, edge_vec &minimal_hitting_sets) {

    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::queue<int> available_workers;
        int num_processes;
        MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

        // Initialize worker pool
        for (int i = 1; i < num_processes; ++i) available_workers.push(i);

        // Convert edges to vectors
        std::vector<int> xv_vec(xv.size()), y_vec(y.size()), x_vec(x.size()), yv_vec(yv.size());
        for (size_t i = 0; i < xv.size(); ++i) xv_vec[i] = xv[i];
        for (size_t i = 0; i < y.size(); ++i) y_vec[i] = y[i];
        for (size_t i = 0; i < x.size(); ++i) x_vec[i] = x[i];
        for (size_t i = 0; i < yv.size(); ++i) yv_vec[i] = yv[i];

		// Serialize minimal_hitting_sets into a flat vector
		std::vector<int> serialized_hitting_sets;
		for (const auto &set : minimal_hitting_sets) {
			for (size_t i = 0; i < set.size(); ++i) {
				serialized_hitting_sets.push_back(set[i] ? 1 : 0);
			}
		}
		// Prepare size metadata
		int hitting_sets_size = serialized_hitting_sets.size();	

        // Assign first task to available workers
        int worker_rank1 = available_workers.front();
        available_workers.pop();

        int worker_rank2 = available_workers.front();
        available_workers.pop();

        // std::cout << "Master (rank 0) assigning task to worker " << worker_rank1 << " for xv and y.\n";
        // std::cout << "Master (rank 0) assigning task to worker " << worker_rank2 << " for x and yv.\n";

        // Create arrays to store MPI_Request objects for the non-blocking sends
		MPI_Request requests[10];
		int req_count = 0;

		// Use MPI_Isend to send tasks and recursion levels for worker_rank1
		MPI_Isend(xv_vec.data(), xv_vec.size(), MPI_INT, worker_rank1, 1, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(y_vec.data(), y_vec.size(), MPI_INT, worker_rank1, 2, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(&hitting_sets_size, 1, MPI_INT, worker_rank1, 3, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, worker_rank1, 4, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(&r, 1, MPI_UNSIGNED_LONG, worker_rank1, 5, MPI_COMM_WORLD, &requests[req_count++]);

		// Use MPI_Isend to send tasks and recursion levels for worker_rank2
		MPI_Isend(x_vec.data(), x_vec.size(), MPI_INT, worker_rank2, 1, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(yv_vec.data(), yv_vec.size(), MPI_INT, worker_rank2, 2, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(&hitting_sets_size, 1, MPI_INT, worker_rank2, 3, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, worker_rank2, 4, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Isend(&r, 1, MPI_UNSIGNED_LONG, worker_rank2, 5, MPI_COMM_WORLD, &requests[req_count++]);

		// Wait for all sends to complete before reusing or modifying the buffers
		MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

		// Dynamic array to hold tasks
		std::vector<WorkData> tasks;

		// Handle results dynamically
		while (!(tasks.empty() && (available_workers.size() == num_processes - 1))) {
			// std::cout << " Available Workers : " << available_workers.size() << "\n";
			// std::cout << " Number of Tasks   : " << tasks.size()<< "\n";

			if (!tasks.empty() && !available_workers.empty()) {
				int worker_rank = available_workers.front();
				available_workers.pop();

				WorkData task = tasks.back();
				tasks.pop_back();

				// Prepare MPI_Request objects for non-blocking sends
				MPI_Request send_requests[5];
				int req_count = 0;

				// Use MPI_Isend to send task data to the worker
				MPI_Isend(task.x_vec.data(), task.x_vec.size(), MPI_INT, worker_rank, 1, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(task.y_vec.data(), task.y_vec.size(), MPI_INT, worker_rank, 2, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&task.hitting_sets_size, 1, MPI_INT, worker_rank, 3, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(task.serialized_hitting_sets.data(), task.hitting_sets_size, MPI_INT, worker_rank, 4, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&task.r, 1, MPI_UNSIGNED_LONG, worker_rank, 5, MPI_COMM_WORLD, &send_requests[req_count++]);

				// Wait for all sends to complete
				MPI_Waitall(req_count, send_requests, MPI_STATUSES_IGNORE);

				// std::cout << "Master sent task to worker " << worker_rank << "\n";
			}

			// Prepare MPI_Request objects for non-blocking receives
			std::vector<MPI_Request> recv_requests(5);
			std::vector<int> recv_x_vec(xv.size()), recv_y_vec(yv.size());
			std::vector<int> serialized_hitting_sets;
			edge::size_type recv_r;
			int hitting_sets_size;

			// Probe to find the source and tag of the incoming message
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			// Get the source of the message
			int source = status.MPI_SOURCE;

			// Use MPI_Irecv for non-blocking receives
			MPI_Irecv(recv_x_vec.data(), recv_x_vec.size(), MPI_INT, source, 1, MPI_COMM_WORLD, &recv_requests[0]);
			MPI_Irecv(recv_y_vec.data(), recv_y_vec.size(), MPI_INT, source, 2, MPI_COMM_WORLD, &recv_requests[1]);
			MPI_Irecv(&hitting_sets_size, 1, MPI_INT, source, 3, MPI_COMM_WORLD, &recv_requests[2]);

			// Wait for the hitting_sets_size to be received before allocating the buffer
			MPI_Wait(&recv_requests[2], MPI_STATUS_IGNORE);
			serialized_hitting_sets.resize(hitting_sets_size);

			MPI_Irecv(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, source, 4, MPI_COMM_WORLD, &recv_requests[3]);
			MPI_Irecv(&recv_r, 1, MPI_UNSIGNED_LONG, source, 5, MPI_COMM_WORLD, &recv_requests[4]);

			// Wait for all the receives to complete
			MPI_Waitall(5, recv_requests.data(), MPI_STATUSES_IGNORE);

			// Check if the worker is signaling that it's free
			if (recv_r == 0) {
				// Worker is free, add it back to the pool of available workers
				available_workers.push(source);
				// std::cout << "Master: Worker " << source << " is free (r = 0).\n";
			} else {
				// Store received data into the task array for further processing
				WorkData new_task = {recv_x_vec, recv_y_vec, hitting_sets_size, serialized_hitting_sets, recv_r};
				tasks.push_back(new_task);
				// std::cout << "Master: Received task result from worker " << source << " with r = " << recv_r << ".\n";
			}
        }

		if (tasks.empty() && (available_workers.size() == num_processes - 1)){
			// Send termination signals to all workers
			std::vector<int> terminate_signal(2, -1); // Termination signal vector (-1 indicates termination)
			std::vector<MPI_Request> termination_requests((num_processes - 1) * 5); // Adjust for 5 sends per worker
			int req_count = 0;

			// Loop through all worker processes (excluding the master)
			for (int i = 1; i < num_processes; ++i) {
				// Send termination signal in the same structure as task assignment
				MPI_Isend(terminate_signal.data(), terminate_signal.size(), MPI_INT, i, 1, MPI_COMM_WORLD, &termination_requests[req_count++]);
				MPI_Isend(terminate_signal.data(), terminate_signal.size(), MPI_INT, i, 2, MPI_COMM_WORLD, &termination_requests[req_count++]);
				MPI_Isend(&terminate_signal[0], 1, MPI_INT, i, 3, MPI_COMM_WORLD, &termination_requests[req_count++]); // Simulating hitting_sets_size
				MPI_Isend(terminate_signal.data(), terminate_signal.size(), MPI_INT, i, 4, MPI_COMM_WORLD, &termination_requests[req_count++]);
				MPI_Isend(&terminate_signal[0], 1, MPI_UNSIGNED_LONG, i, 5, MPI_COMM_WORLD, &termination_requests[req_count++]); // Simulating r
			}

			// Wait for all termination signals to be sent
			MPI_Waitall(req_count, termination_requests.data(), MPI_STATUSES_IGNORE);

			// std::cout << " >>> Sent Termination signals to all workers.\n";
		}			

    } else {
        worker();  // Non-master processes become workers
    }
	MPI_Finalize();
	return;
}

void Hypergraph::worker() {
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::cout << "Worker " << rank << " started and waiting for tasks from master (rank 0).\n";

    bool extend = false;
    std::vector<int> a_vec(m_num_vertices);
    std::vector<int> b_vec(m_num_vertices);
    edge::size_type r;

    int hitting_sets_size;
    edge_vec minimal_hitting_sets;
    std::vector<int> serialized_hitting_sets;
    
    edge a;
    edge b;
    bool terminate;

    while (true) {

        if (!extend) {
            terminate = true;

            MPI_Request recv_requests[5];

            // Non-blocking receive for a_vec
            MPI_Irecv(a_vec.data(), a_vec.size(), MPI_INT, 0, 1, MPI_COMM_WORLD, &recv_requests[0]);

            // Non-blocking receive for b_vec
            MPI_Irecv(b_vec.data(), b_vec.size(), MPI_INT, 0, 2, MPI_COMM_WORLD, &recv_requests[1]);

            // Non-blocking receive for hitting_sets_size
            MPI_Irecv(&hitting_sets_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &recv_requests[2]);

            // Non-blocking receive for serialized_hitting_sets (placeholder size initially)
            MPI_Irecv(serialized_hitting_sets.data(), serialized_hitting_sets.size(), MPI_INT, 0, 4, MPI_COMM_WORLD, &recv_requests[3]);

            // Non-blocking receive for r
            MPI_Irecv(&r, 1, MPI_UNSIGNED_LONG, 0, 5, MPI_COMM_WORLD, &recv_requests[4]);

            // Wait for a_vec to check for termination signal
            MPI_Wait(&recv_requests[0], &status);

            // Check for termination signal in a_vec
            if (!a_vec.empty() && a_vec[0] == -1) {
                // std::cout << "#### Worker " << rank << " received Termination signal.\n";
                break; // Exit the worker loop immediately
            }

            // Wait for hitting_sets_size
            MPI_Wait(&recv_requests[2], &status);

            // Resize serialized_hitting_sets based on hitting_sets_size
            serialized_hitting_sets.resize(hitting_sets_size);

            // Wait for all other receives to complete
            MPI_Waitall(4, &recv_requests[1], MPI_STATUSES_IGNORE);

            // Deserialize minimal_hitting_sets
            minimal_hitting_sets.clear();
            for (size_t i = 0; i < hitting_sets_size; i += m_num_vertices) {
                edge set(m_num_vertices);
                for (size_t j = 0; j < m_num_vertices; ++j) {
                    set[j] = serialized_hitting_sets[i + j];
                }
                minimal_hitting_sets.push_back(set);
            }


			// Convert vectors to edges safely
			a.resize(a_vec.size());
			a.reset();
			b.resize(b_vec.size());
			b.reset();

			// Convert vectors to edges
			for (size_t i = 0; i < a_vec.size(); ++i) a[i] = a_vec[i];
			for (size_t i = 0; i < b_vec.size(); ++i) b[i] = b_vec[i];
        }


		// // Debug: Print task received by worker
        // std::cout << "Worker " << rank << " received task with vectors:\n";
        // std::cout << "a_vec: ";
        // for (auto v : a_vec) std::cout << v << " ";
        // std::cout << "\nb_vec: ";
        // for (auto v : b_vec) std::cout << v << " ";
        // std::cout << std::endl;
		// std::cout << "\nr: ";
		// std::cout << r <<std::endl;

        // Check extendability
		if(terminate){
        	extend = extendable(a, b);
		}

        if (extend) {
            if (r == m_num_vertices) {
                minimal_hitting_sets.push_back(a);  // Correct edge assignment
				// std::cout << "\nWorker " << rank << " found set --  ";
				// std::cout << "\nWorker " << rank << " terminating with minimal hitting set: ";
				// for (size_t i = 0; i < a.size(); ++i) {
				// 	std::cout << a[i] << " ";
				// }
				// std::cout << "\n";

                if (m_configuration.collect_hitting_set_statistics) {
                    auto now = Clock::now();
                    m_hitting_set_stats.add_record({ edge_to_string(a), ns_string(m_hitting_set_timestamp, now) });
                    m_hitting_set_timestamp = Clock::now();
                }
				terminate = false;
				extend = false;

				MPI_Request send_requests[5];
				int req_count = 0;
				r = 0; // Indicate worker is free

				// Non-blocking send to master
				MPI_Isend(a_vec.data(), a_vec.size(), MPI_INT, 0, 1, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(b_vec.data(), b_vec.size(), MPI_INT, 0, 2, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&hitting_sets_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, 0, 4, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&r, 1, MPI_UNSIGNED_LONG, 0, 5, MPI_COMM_WORLD, &send_requests[req_count++]);

				// Wait for all sends to complete
				MPI_Waitall(req_count, send_requests, MPI_STATUSES_IGNORE);
            }else {
				// Proper assignment with safer Boost functions
				edge xv = a, yv = b;
				xv.set(r); 
				yv.set(r);
				r++;

				std::vector<int> x_vec(a.size()), yv_vec(yv.size());
				for (size_t i = 0; i < a.size(); ++i) x_vec[i] = a[i];
				for (size_t i = 0; i < yv.size(); ++i) yv_vec[i] = yv[i];

				// std::cout << "xv: ";
				// for (auto v : x_vec) std::cout << v << " ";
				// std::cout << "\nyv: ";
				// for (auto v : yv_vec) std::cout << v << " ";
				// std::cout << std::endl;

				// Prepare MPI_Request objects for non-blocking sends
				MPI_Request send_requests[5];
				int req_count = 0;

				// Use MPI_Isend to send data back to the master
				MPI_Isend(x_vec.data(), x_vec.size(), MPI_INT, 0, 1, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(yv_vec.data(), yv_vec.size(), MPI_INT, 0, 2, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&hitting_sets_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, 0, 4, MPI_COMM_WORLD, &send_requests[req_count++]);
				MPI_Isend(&r, 1, MPI_UNSIGNED_LONG, 0, 5, MPI_COMM_WORLD, &send_requests[req_count++]);

				// Wait for all sends to complete
				MPI_Waitall(req_count, send_requests, MPI_STATUSES_IGNORE);

				// std::cout << "Worker " << rank << " sent task successfully.\n";

				// Continue processing with a = xv
				a = xv;

				// Convert edge `a` back to a vector for debugging
				std::vector<int> a_vec(a.size());
				for (size_t i = 0; i < a.size(); ++i) a_vec[i] = a[i];

				// Debug output for `a`
				// std::cout << "\n a = xv: ";
				// for (const auto& v : a_vec) std::cout << v << " ";
				// std::cout << std::endl;
			}
        }else{
			MPI_Request send_requests[5];
			int req_count = 0;
			r = 0; // Indicate worker is free

			// Non-blocking send to master
			MPI_Isend(a_vec.data(), a_vec.size(), MPI_INT, 0, 1, MPI_COMM_WORLD, &send_requests[req_count++]);
			MPI_Isend(b_vec.data(), b_vec.size(), MPI_INT, 0, 2, MPI_COMM_WORLD, &send_requests[req_count++]);
			MPI_Isend(&hitting_sets_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &send_requests[req_count++]);
			MPI_Isend(serialized_hitting_sets.data(), hitting_sets_size, MPI_INT, 0, 4, MPI_COMM_WORLD, &send_requests[req_count++]);
			MPI_Isend(&r, 1, MPI_UNSIGNED_LONG, 0, 5, MPI_COMM_WORLD, &send_requests[req_count++]);

			// Wait for all sends to complete
			MPI_Waitall(req_count, send_requests, MPI_STATUSES_IGNORE);

			// std::cout << "Worker " << rank << " is free (r = 0) and reported back to master.\n";
		}
    }
	return;
}



void Hypergraph::minimize() {
	edge_set new_edges;
	for (edge old_edge : m_edges) {
		bool insert = true;
		edge_vec to_delete;
		for (edge new_edge : new_edges) {
			if (new_edge.is_subset_of(old_edge)) {
				insert = false;
				break;
			}
			if (old_edge.is_subset_of(new_edge)) { // edge is subset of other one
				to_delete.push_back(new_edge);
			}
		}
		if (insert) {
			for (auto e : to_delete) new_edges.erase(e);
			new_edges.insert(old_edge);
		}
	}
	m_edges = edge_vec(new_edges.begin(), new_edges.end());
}

void Hypergraph::permute(permutation p) {
	if (m_edges.empty()) return;
	if (p.size() != m_num_vertices) {
		std::cerr << "Cannot apply permutation of length " << p.size() << " to graph with " << m_num_vertices << " vertices!" << std::endl;
		return;
	}
	edge_vec new_edges;
	for (edge e : m_edges) {
		edge new_edge(m_num_vertices);
		for (auto i = e.find_first(); i != edge::npos; i = e.find_next(i)) new_edge[p[i]] = 1;
		new_edges.push_back(new_edge);
	}
	m_edges = new_edges;
}

edge_vec Hypergraph::brute_force_mhs() {
	// brute force all minimal hitting sets
	// this works a bit similar to the apriori algorithm
	int check_time_and_memory_counter = 0;
	auto exit_timestamp = Clock::now();
	if (m_configuration.collect_hitting_set_statistics) m_hitting_set_timestamp = Clock::now();
	edge_vec minimal_hitting_sets;
	edge_vec incomplete_hitting_sets;
	for (edge::size_type i = 0; i < m_num_vertices; ++i) {
		edge e(m_num_vertices);
		e[i] = 1;
		if (is_hitting_set(e)) {
			minimal_hitting_sets.push_back(e);
		}
		else incomplete_hitting_sets.push_back(e);
	}
	if (incomplete_hitting_sets.empty()) return minimal_hitting_sets;
	for (int set_size = 1; set_size < m_num_vertices; ++set_size) {
		edge_vec::size_type current_level_cutoff = minimal_hitting_sets.size();
		edge_vec new_incomplete_hitting_sets;
		for (edge_vec::size_type i_first_set = 0; i_first_set < incomplete_hitting_sets.size(); ++i_first_set) {
			for (edge_vec::size_type i_second_set = i_first_set + 1; i_second_set < incomplete_hitting_sets.size(); ++i_second_set) {
				edge first_set_minus_last_vertex = incomplete_hitting_sets[i_first_set];
				first_set_minus_last_vertex[last_vertex(first_set_minus_last_vertex)] = 0;
				if (first_set_minus_last_vertex.is_subset_of(incomplete_hitting_sets[i_second_set])) {
					edge candidate = incomplete_hitting_sets[i_first_set] | incomplete_hitting_sets[i_second_set];
					bool add_candidate = true;
					for (edge_vec::size_type i = 0; i < current_level_cutoff; ++i) {
						if (minimal_hitting_sets[i].is_subset_of(candidate)) {
							add_candidate = false;
							break;
						}
					}
					if (add_candidate) {
						if (is_hitting_set(candidate)) {
							minimal_hitting_sets.push_back(candidate);
							if (minimal_hitting_sets.size() == 81) {
								std::cout << "";
							};
							if (m_configuration.collect_hitting_set_statistics) {
								auto now = Clock::now();
								m_hitting_set_stats.add_record({ edge_to_string(candidate), ns_string(m_hitting_set_timestamp, now) });
								m_hitting_set_timestamp = Clock::now();
							}
						}
						else new_incomplete_hitting_sets.push_back(candidate);
					}
				}
				else break;
			}
			check_time_and_memory_counter++;
			if (check_time_and_memory_counter > 1000) {
				check_time_and_memory_counter = 0;
				if (std::chrono::duration_cast<std::chrono::minutes>(Clock::now() - exit_timestamp).count() >= 60 * 12) {
					std::cerr << "Aborting due to time constraints, took " << time_string(exit_timestamp, Clock::now()) << "." << std::endl;
					return edge_vec();
				}				
			}
		}
		if (new_incomplete_hitting_sets.empty()) return minimal_hitting_sets;
		incomplete_hitting_sets.swap(new_incomplete_hitting_sets);
	}
	return minimal_hitting_sets;
}

int Hypergraph::maximum_iteration_count(std::vector<edge_vec> s) {
	if (s.empty()) return 0;
	int maximum_iteration_count = 1;
	for (edge_vec sx : s)  maximum_iteration_count *= sx.size();
	return maximum_iteration_count;
}

int Hypergraph::summed_sx_sizes(std::vector<edge_vec> s) {
	int summed_sx_sizes = 0;
	for (edge_vec sx : s) summed_sx_sizes += sx.size();
	return summed_sx_sizes;
}

int Hypergraph::total_number_of_vertices_in_s(std::vector<edge_vec> s) {
	int total_number_of_vertices_in_s = 0;
	for (edge_vec sx: s) for (edge e : sx) total_number_of_vertices_in_s += e.count();
	return total_number_of_vertices_in_s;
}

int Hypergraph::total_number_of_vertices_in_t(edge_vec t) {
	int total_number_of_vertices_in_t = 0;
	for (edge e : t) total_number_of_vertices_in_t += e.count();
	return total_number_of_vertices_in_t;
}

void Hypergraph::save_statistics() {
	if (m_configuration.collect_hitting_set_statistics) {
		fs::path hitting_set_statistics_path = m_configuration.statistics_directory;
		hitting_set_statistics_path /= (m_configuration.name + "_" + m_configuration.implementation + "_hitting_set_statistics.csv");
		m_hitting_set_stats.save(hitting_set_statistics_path.string());
	}
	if (m_configuration.collect_oracle_statistics) {
		fs::path oracle_statistics_path = m_configuration.statistics_directory;
		oracle_statistics_path /= (m_configuration.name + "_" + m_configuration.implementation + "_oracle_statistics.csv");
		m_oracle_stats.save(oracle_statistics_path.string());
	}
}