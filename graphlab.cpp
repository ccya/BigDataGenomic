#include <string>
#include <graphlab.hpp>

// Define the graph
struct reads {
    std::string readname;
    std:string content;
    double score
    reads():score(0.0) { }
    explicit reads(std::string name, std::string cont):readname(name), content(comt),score(0.0){ }
    
    // To make it serializable
    void save(graphlab::oarchive& oarc) const {
      oarc << readname << content << score;
    }
    void load(graphlab::iarchive& iarc) {
      iarc >> readname >> content >> score;
    }
  };
 	// Define the type of graph
  typedef graphlab::distributed_graph<reads, graphlab::empty> graph_type;
  
// Load the graph data
    bool line_parser(graph_type& graph, 
                   const std::string& filename, 
                   const std::string& textline) {
    std::stringstream strm(textline);
    graphlab::vertex_id_type vid;
    std::string readname;
    std::string content;
    // first entry in the line is a vertex ID
    strm >> vid;
    strm >> readname;
    strm >> content;
    // insert this web page
    graph.add_vertex(vid, reads(readname,content));
    // while there are elements in the line, continue to read until we fail
    while(1){
      graphlab::vertex_id_type other_vid;
      strm >> other_vid;
      if (strm.fail()) 
        return false;
      graph.add_edge(vid, other_vid);
    }
    retun true;
  }



  int main(int argc, char** argv) {
    graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
 
    graph_type graph(dc);
    graph.load("graph.txt", line_parser);
    graphlab::mpi_tools::finalize();
  }