# Federated Enumhyp UCC Discovery

A research-oriented Java framework for **federated discovery of Unique Column Combinations (UCCs)** using a **server–client execution model**, **containerized local enumeration**, and **centralized aggregation of derived results**. The project investigates whether minimal UCC discovery can be executed efficiently when raw tuples remain local to participating nodes and only derived UCC candidates are shared.

This repository accompanies the implementation and evaluation described in the project report **“Federated Data Profiling and Unique Column Combination Discovery: Decentralized Execution with Central Coordination”**.

---

## 1. Project Motivation

Unique Column Combinations (UCCs) are minimal sets of attributes that uniquely identify tuples in a relational dataset. UCC discovery is a core task in data profiling, schema analysis, and data quality assessment. In classical settings, UCC discovery is performed centrally on the full dataset. In privacy-sensitive or distributed environments, however, raw tuples may not be transferable to a central coordinator.

This project addresses that constraint by studying a **federated alternative**:

- each client processes its own data partition locally,
- local UCCs are computed through the `enumhyp` pipeline,
- only the resulting column-index combinations are transmitted,
- a server aggregates local outputs into a global result.

The underlying research question is whether this strategy can preserve usefulness while reducing end-to-end runtime on larger datasets.

---

## 2. Research Context

The implementation builds on the `enumhyp` toolchain and the hitting-set based formulation of UCC discovery. In this formulation:

- columns correspond to vertices,
- difference sets between tuples correspond to hyperedges,
- minimal hitting sets correspond to minimal UCCs.

The broader algorithmic foundation comes from the efficient enumeration of minimal hitting sets with decision-tree exploration and extension-oracle pruning. This repository does **not** reimplement the enumeration algorithm from scratch in Java; instead, it orchestrates and integrates the external `enumhyp` workflow inside a federated execution setup.

---

## 3. What This Repository Implements

The repository provides an end-to-end prototype for comparing:

1. **Centralized baseline execution** on the full dataset.
2. **Federated execution** across multiple client partitions.
3. **Server-side aggregation** of local UCC outputs.
4. **Result comparison** between distributed and centralized outputs.

The implementation emphasizes:

- reproducible execution through Docker,
- strict separation between local data and shared derived results,
- modular Java code for orchestration, networking, file handling, and aggregation,
- experimental evaluation of runtime tradeoffs.

---

## 4. System Architecture

The system follows a **decentralized execution with central coordination** design.

### 4.1 Client-side processing

Each federated node receives one partition of a dataset and performs the following steps locally:

1. read a local CSV partition,
2. execute the `generate` stage to construct a UCC graph,
3. execute the `enumerate` stage to derive minimal UCCs,
4. normalize and send only derived UCC results to the server.

No raw tuple values are transmitted.

### 4.2 Server-side coordination

The coordinator:

1. accepts incoming client connections,
2. stores each client's transmitted results separately,
3. aggregates the received local UCC sets,
4. writes the resulting global UCCs and logs to output files.

### 4.3 Aggregation strategies

Two aggregation strategies are considered:

#### Intersection-based aggregation
Returns only UCCs that occur in **all** local result sets.

**Strengths**
- simple,
- fast,
- practical in identical-schema partitioning scenarios.

**Limitations**
- strict,
- may miss valid global UCCs that are not repeated verbatim across all nodes.

#### Union-based MinUnion aggregation
Builds global candidates by taking unions of local UCCs and keeping only subset-minimal combinations.

**Strengths**
- better aligned with the union-of-constraints interpretation.

**Limitations**
- combinatorially expensive,
- unsuitable as a default choice for larger client counts or large local UCC sets.

In the reported experiments, **intersection** is the practical default because it is computationally lightweight and stable for the chosen split setting.

---

## 5. Repository Structure

The visible repository structure indicates the following main components:

```text
.
├── data/
├── serverReceived/
├── src/
│   ├── main/java/at/univie/
│   │   ├── centralized/
│   │   │   └── CentralizedBaseline.java
│   │   ├── enumhyp/
│   │   │   ├── EnumhypRunner.java
│   │   │   └── GenerateEnumerate.java
│   │   └── federated/
│   │       ├── aggregator/
│   │       │   ├── FilePathManager.java
│   │       │   ├── UccAggregator.java
│   │       │   ├── UccFileReader.java
│   │       │   └── UccMinimizer.java
│   │       ├── client/
│   │       │   └── Client.java
│   │       ├── compared/
│   │       │   └── FileComparator.java
│   │       ├── server/
│   │       │   ├── ErrorHandler.java
│   │       │   ├── FileManager.java
│   │       │   ├── MessageHandler.java
│   │       │   └── Server.java
│   │       └── main.java
│   └── test/java/at/univie/federated/
├── Dockerfile
└── pom.xml
```

### 5.1 Directory purpose

- `data/` stores dataset partitions and generated intermediate artifacts.
- `serverReceived/` stores per-client uploads received by the coordinator.
- `centralized/` contains the reference baseline execution logic.
- `enumhyp/` encapsulates execution of the external Docker-based enumeration pipeline.
- `federated/client/` contains the client-side networking and transmission logic.
- `federated/server/` contains the TCP server, message handling, and persistence logic.
- `federated/aggregator/` contains aggregation and minimality reduction utilities.
- `federated/compared/` provides comparison logic between baseline and federated outputs.

---

## 6. Technology Stack

- **Java 11**
- **Maven**
- **Docker**
- **Ubuntu-based enumhyp container**
- **TCP socket communication** for client–server coordination

The `pom.xml` declares a Java 11 Maven build, and the `Dockerfile` builds an Ubuntu image that clones and compiles `enumhyp` from source.

---

## 7. Dockerized Enumhyp Pipeline

The project uses a dedicated Docker image to provide a uniform runtime environment for local UCC discovery. The container:

- starts from `ubuntu:22.04`,
- installs build tools and Boost dependencies,
- clones the `enumhyp` implementation,
- compiles it via CMake,
- exposes the compiled binary as the container entrypoint.

This design removes host-side dependency drift and ensures that every node executes the same enumeration environment.

---

## 8. Build Instructions

### 8.1 Prerequisites

Ensure the following are available on your machine:

- Java 11
- Maven 3.x
- Docker
- network access for Docker image creation and dependency retrieval

### 8.2 Build the Java project

```bash
mvn clean package
```

### 8.3 Build the Docker image for enumhyp

```bash
docker build -t enumhyp-ubuntu .
```

---

## 9. Expected Execution Workflow

Because this repository is structured as a research prototype rather than a polished CLI application, execution is best understood as a staged workflow.

### Step 1 — Prepare data partitions
Place dataset partitions in the `data/` directory using a consistent naming scheme such as:

```text
flight_p1.csv
flight_p2.csv
...
flight_pn.csv
```

### Step 2 — Run local UCC discovery
Each client invokes the Dockerized `generate` and `enumerate` stages on its local partition and produces a local text file containing minimal UCCs represented as column-index sets.

### Step 3 — Start the server
Launch the coordinator so it can accept client connections and store incoming UCC streams in `serverReceived/`.

### Step 4 — Launch clients
Each client connects to the server and transmits its derived local UCC output line by line.

### Step 5 — Aggregate local UCCs
After all client uploads are complete, the server-side aggregation step derives the final global UCC list.

### Step 6 — Compare with centralized baseline
Run the centralized pipeline on the full dataset and compare outputs and runtimes.

---

## 10. Input and Output Conventions

### Input

Typical inputs include:

- CSV partitions stored locally at each client,
- shared-schema datasets split from one original table,
- optional baseline full-dataset CSV file.

### Intermediate artifacts

The execution pipeline produces artifacts such as:

- generated UCC graph files,
- local UCC text files,
- per-client received files on the server.

### Output

Representative output files include:

- per-client UCC files,
- `serverpart_*.txt` files in `serverReceived/`,
- aggregated global UCC output,
- comparison artifacts between baseline and federated results.

---

## 11. Experimental Findings

The project report evaluates the system on three datasets: **abalone**, **fdReduced**, and **flight**.

The main findings are:

- For **small datasets**, federated execution may be slightly slower because container startup, communication, and coordination overhead dominate the workload.
- For **larger datasets**, federated execution substantially reduces end-to-end runtime because local enumeration is distributed across multiple nodes.
- **Intersection-based aggregation** remains negligible in cost and is operationally practical.
- **Union-based MinUnion aggregation** becomes a scalability bottleneck due to combinatorial growth.

Reported runtime observations from the project report include:

- `fdReduced (n = 10)`: approximately **2.70 min → 0.52 min** from centralized to distributed execution.
- `flight (n = 20)`: approximately **40.87 min → 4.01 min** from centralized to distributed execution.
- Union aggregation on `flight (n = 20)` reaches an outlier of approximately **16.99 min**.

These results indicate that the proposed architecture is most effective when local enumeration dominates the total runtime and the server-side aggregation remains lightweight.

---

## 12. Design Strengths

- Preserves local data ownership by avoiding raw tuple transfer.
- Uses Docker for reproducibility across nodes.
- Separates concerns cleanly across baseline, client, server, and aggregation modules.
- Demonstrates clear runtime gains on larger datasets.
- Provides a useful experimental bridge between centralized data profiling and federated execution.

---

## 13. Current Limitations

This repository is a **research prototype**, not a production-ready federated data profiling platform.

Known limitations include:

- aggregation correctness is approximation-sensitive, especially for intersection,
- union-based aggregation does not scale well,
- execution appears class-driven rather than exposed through a unified command-line interface,
- deployment orchestration across multiple hosts is manual,
- privacy guarantees are architectural rather than formally cryptographic.

---

## 14. Reproducibility Notes

For reproducible experiments:

1. use the same Java and Docker versions across machines,
2. keep dataset schemas identical across partitions,
3. normalize local UCC outputs before aggregation,
4. record client counts, dataset sizes, and runtime measurements consistently,
5. compare every federated run against a centralized reference result.

---

## 15. Suggested Future Extensions

This repository is a solid foundation for several research directions:

- stronger aggregation methods with correctness guarantees,
- privacy-enhanced result sharing,
- orchestration via Docker Compose or Kubernetes,
- automated benchmarking and result logging,
- support for heterogeneous nodes and failure handling,
- richer evaluation on real-world distributed data profiling workloads.

---

## 16. Citation-Style Project Summary

If you reference this repository in academic or project documentation, the following summary is appropriate:

> A Java-based federated prototype for Unique Column Combination discovery that combines Dockerized `enumhyp` execution at local nodes with server-side aggregation of derived UCC results, enabling empirical comparison between centralized and distributed execution strategies.

---

## 17. Acknowledgement of Foundations

This project is based on the hitting-set interpretation of UCC discovery and integrates the external `enumhyp` implementation into a federated orchestration framework. The implementation and evaluation are closely aligned with the accompanying University of Vienna Praktikum project work.

---

## 18. References

1. T. Bläsius, T. Friedrich, J. Lischeid, K. Meeks, M. Schirneck. *Efficiently Enumerating Hitting Sets of Hypergraphs Arising in Data Profiling*. Journal of Computer and System Sciences, 2022.
2. Hasso-Plattner-Institut. *Enumeration in Data Profiling*.
3. `enumhyp` implementation by Julius Lischeid / related project sources.

---
