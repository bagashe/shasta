#include "MetaMarkerGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>


void MetaMarkerGraph::createEdges()
{
    using Graph = MetaMarkerGraph;
    Graph& graph = *this;



    // Construct the sequence of vertices encountered by each oriented read.
    std::map<OrientedReadId, vector<vertex_descriptor> > pseudoPaths;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const vector< pair<OrientedReadId, uint64_t> >& orientedReads =
            graph[v].orientedReads;

        // Loop over oriented reads of this vertex and the corresponding
        // metaOrdinals.
        for(const auto& p: orientedReads) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t metaOrdinal = p.second;

            // Access the pseudo-path for this oriented read.
            vector<vertex_descriptor>& pseudoPath = pseudoPaths[orientedReadId];

            // Make sure we have a slot for this metaOrdinal.
            if(pseudoPath.size() <= metaOrdinal) {
                pseudoPath.resize(metaOrdinal+1, null_vertex());
            }

            // Now we can store it.
            pseudoPath[metaOrdinal] = v;
        }
    }


    // Check that the pseudo-path don't have any missing vertices.
    for(const auto& p: pseudoPaths) {
        const vector<vertex_descriptor>& pseudoPath = p.second;
        for(const vertex_descriptor v: pseudoPath) {
            SHASTA_ASSERT(v != null_vertex());
        }
    }




    // Now we can create the edges by looping over the pseudo-path
    // of each orientedRead.
    for(const auto& p: pseudoPaths) {
        const OrientedReadId orientedReadId = p.first;
        const vector<vertex_descriptor>& pseudoPath = p.second;

        // Loop over successive vertices in the pseudo-path of this oriented read.
        for(uint64_t metaOrdinal1=1; metaOrdinal1<pseudoPath.size(); metaOrdinal1++) {
            const uint64_t metaOrdinal0 = metaOrdinal1 - 1;
            const vertex_descriptor v0 = pseudoPath[metaOrdinal0];
            const vertex_descriptor v1 = pseudoPath[metaOrdinal1];

            // Access the edge between these two vertices, creating it if necessary.
            bool edgeExists = false;;
            edge_descriptor e;
            tie(e, edgeExists) = edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = add_edge(v0, v1, graph);
            }

            // Store this metaOrdinal in the edge.
            graph[e].orientedReads.push_back(make_pair(orientedReadId, metaOrdinal0));
        }
    }
}



void MetaMarkerGraph::writeGraphviz(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream graphOut(fileName);
    graphOut << "digraph MetaMarkerGraph {\n";
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        graphOut << graph[v0].vertexId << "->" << graph[v1].vertexId <<
            " [penwidth=" << int(0.3* double(graph[e].orientedReads.size())) << "]"
            << ";\n";
    }
    graphOut << "}\n";
}


void MetaMarkerGraph::writeGfa(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";


    // Write a segment record for each vertex.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        gfa <<
            "S\t" <<
            graph[v].segmentId << "\t" <<
            "*\t" <<
            "LN:i:" << graph[v].markerCount <<
            "\n";
    }

    // Write link records.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        BGL_FORALL_OUTEDGES(v0, e01, graph, Graph) {
            const vertex_descriptor v1 = target(e01, graph);
            gfa <<
                "L\t" <<
                graph[v0].vertexId << "\t" <<
                "+\t" <<
                graph[v1].vertexId << "\t" <<
                "+\t" <<
                "*\n";
        }
    }
}



void MetaMarkerGraph::writeVerticesCsv(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId,Segment id,Marker count,Coverage,Segment id and coverage\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const MetaMarkerGraphVertex& vertex = graph[v];
        csv << vertex.vertexId << ",";
        csv << vertex.segmentId << ",";
        csv << vertex.markerCount << ",";
        csv << vertex.orientedReads.size() << ",";
        csv << vertex.segmentId << "/";
        csv << vertex.orientedReads.size() << "\n";
    }

}



void MetaMarkerGraph::writeEdgesCsv(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0,VertexId1,Coverage\n";

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << graph[v0].vertexId << ",";
        csv << graph[v1].vertexId << ",";
        csv << graph[e].orientedReads.size() << "\n";
    }

}

