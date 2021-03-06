#ifndef SHASTA_WRITE_GRAPH_HPP
#define SHASTA_WRITE_GRAPH_HPP

// Code to write a Boost graph directly, without using Graphviz rendering.

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "iostream.hpp"
#include <limits>
#include <map>
#include "string.hpp"

namespace shasta {

    namespace WriteGraph {

        class VertexAttributes {
        public:
            double radius = 1.;
            string id;
            string color = "black";
            string tooltip;
            string url;
        };

        class EdgeAttributes {
        public:
            double thickness = 1.;
            string id;
            string color = "black";
            string tooltip;
            string url;
        };

        template<class Graph> void writeSvg(
            const Graph&,
            const string& svgId,
            uint64_t width,
            uint64_t height,
            const std::map<typename Graph::vertex_descriptor, VertexAttributes>&,
            const std::map<typename Graph::edge_descriptor, EdgeAttributes>&,
            ostream&);
    }
}



template<class Graph> void shasta::WriteGraph::writeSvg(
    const Graph& graph,
    const string& svgId,
    uint64_t width,
    uint64_t height,
    const std::map<typename Graph::vertex_descriptor, VertexAttributes>& vertexAttributes,
    const std::map<typename Graph::edge_descriptor, EdgeAttributes>& edgeAttributes,
    ostream& svg)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    // using edge_descriptor = typename Graph::edge_descriptor;

    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const auto& position = graph[v].position;
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if(it != vertexAttributes.end()) {
            attributes = it->second;
        }
        const double radius = attributes.radius;

        // Update the view box to include this vertex.
        xMin = min(xMin, position[0] - radius);
        xMax = max(xMax, position[0] + radius);
        yMin = min(yMin, position[1] - radius);
        yMax = max(yMax, position[1] + radius);
    }



    // Begin the svg.
    svg << "<svg id='" << svgId << "' width='" << width << "' height='" << height <<
        "' viewbox='" << xMin << " " << yMin << " " << xMax-xMin << " " << yMax-yMin <<
        "'>\n";



    // Write the edges first, so they don't cover the vertices.
    svg << "<g id='" << svgId << "-edges'>\n";
    BGL_FORALL_EDGES_T(e, graph, Graph) {

        // Get the attributes for this vertex.
        EdgeAttributes attributes;
        auto it = edgeAttributes.find(e);
        if(it != edgeAttributes.end()) {
            attributes = it->second;
        }

        // Get vertex positions.
        const vertex_descriptor v1 = source(e, graph);
        const vertex_descriptor v2 = target(e, graph);
        const auto& position1 = graph[v1].position;
        const auto& position2 = graph[v2].position;

        svg << "<line x1='" << position1[0] << "' y1='" << position1[1] <<
            "' x2='" << position2[0] << "' y2='" << position2[1];

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << "' stroke='" << attributes.color <<
            "' stroke-width='" << attributes.thickness <<
            "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</line>\n";
    }
    svg << "</g>\n";



    // Write the vertices.
    svg << "<g id='" << svgId << "-vertices' stroke='none'>\n";
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if(it != vertexAttributes.end()) {
            attributes = it->second;
        }
        const auto& position = graph[v].position;

        if(not attributes.url.empty()) {
            svg << "<a href='" << attributes.url << "'>";
        }

        svg << "<circle cx='" << position[0] << "' cy='" << position[1] <<
            "' r='" << attributes.radius << "'";

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << " fill='" << attributes.color << "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</circle>";

        if(not attributes.url.empty()) {
            svg << "</a>";
        }
        svg << "\n";
    }
    svg << "</g>\n";



    // End the svg.
    svg << "</svg>\n";
}


#endif
