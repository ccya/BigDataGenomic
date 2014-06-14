#include <graphlab.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


struct read_data
{	
	int length;
	std::string content;
	int offset;
	double score;
	int end;
    bool valid;
    read_data(): length(0), content(""), offset(0), score(0), end(0), valid(false){ }
	explicit read_data(int length, std::string content, int offset, double score, int end, bool valid) : length(length), content(content), offset(offset), score(score), end(end), valid(valid) { } 
    void save(graphlab::oarchive& oarc) const {
      oarc << length << content << offset << score << end << valid;
    }
    void load(graphlab::iarchive& iarc) {
      iarc >> length >> content >> offset >> score >> end >> valid;
    }
};

// The type of graph used in this program
typedef graphlab::distributed_graph<read_data, graphlab::empty> graph_type;


bool load_graph(const std::string& filename, graph_type& graph) {

	std::ifstream fin(filename.c_str());
    std::cout << "Load " << filename << " now ... " << std::endl;
	if(!fin.good()) return false;

	std::string textline;
	if (fin.is_open()) {
        std::cout << "File is open! " << std::endl;
		graphlab::vertex_id_type vid;
		while (std::getline(fin, textline)) {
            std::cout << "Read a line. " << std::endl;
			std::stringstream strm(textline);

    		int length;
    		std::string content;
    		int offset;
    		double score;
    		//double dist;

    		// first entry in the line is a vertex ID
    		strm >> vid;
    		strm >> length;
    		strm >> content;
    		strm >> offset;
    		strm >> score;
    		int end = offset + length;
            std::cout << "textline: " << textline << std::endl;
    		// insert this read
   			if(graph.add_vertex(vid, read_data(length, content, offset, score, end, true)))
                std::cout << "add_vertex success!" << std::endl;

            // if(graph.contains_vertex(vid)) {std::cout << "gaph contains the vid!!!!!" << std::endl;}
            // if(graph.num_vertices() > 0) {std::cout << last_vid << " not empty" << std::endl;}
         
            // std::cout << "xxxxx " << graph.vertex(vid).id() << std::endl;
            std::cout << "vertex id: " << vid << ", score: " << score << std::endl;

   			// while there are elements in the line, continue to read until we fail
   			while(1){
   				graphlab::vertex_id_type other_vid;
    			strm >> other_vid;
                std::cout << "xxxxxxxxx: " << other_vid << std::endl;
                if (strm.fail())
                    break;
                std::cout << "edge: (" << vid << ", " << other_vid << ")" << std::endl;
    			graph.add_edge(vid, other_vid);
    		}
		}
	}

	return true;
}

// Get the other vertex in the edge.
graph_type::vertex_type get_other_vertex(const graph_type::edge_type& edge, 
                                        const graph_type::vertex_type& vertex) {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
}


// Exempt reads that are leaves but are not the last read
class exempt_reads_program : public graphlab::ivertex_program<graph_type, graphlab::empty, graphlab::vertex_id_type>,
                             public graphlab::IS_POD_TYPE{

    graphlab::vertex_id_type last_id;
public:
    void init(icontext_type& context, const vertex_type& vertex, const graphlab::vertex_id_type& msg) { 
        last_id = msg;
    }

    edge_dir_type gather_edges (icontext_type& context, const vertex_type& vertex) const {
        return graphlab:: NO_EDGES;
    }

    void apply(icontext_type& context, vertex_type& vertex, const graphlab::empty& empty){
        // if a vertex do not have a child and is not the last vertex, then tag it as an invalid vertex.
        if (vertex.num_out_edges() == 0 && vertex.id() != last_id) {
            vertex.data().valid = false;
        }
    }
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return graphlab::NO_EDGES; 
    } // end of scatter_edge
};


// gather type
struct id_and_score
{
    std::vector<graphlab::vertex_id_type> ids;
    std::vector<double> scores;
    id_and_score(): ids(),scores() { }
    id_and_score(graphlab::vertex_id_type id, double score): ids(), scores() {
        ids.push_back(id);
        scores.push_back(score);
    }
    id_and_score& operator += (const id_and_score& other) {
        for (size_t i = 0; i < other.ids.size(); ++i) {
            ids.push_back(other.ids[i]);
            scores.push_back(other.scores[i]);
        }
        return *this;
    }

    void save(graphlab::oarchive& oarc) const {
        size_t num = ids.size();
        oarc << num;
        for (size_t a = 0; a < num; ++a) 
            oarc << ids[a] << scores[a];
    }

    void load(graphlab::iarchive& iarc) {
        ids.clear();
        scores.clear();
        size_t num = 0;
        iarc >> num;
        for (size_t a = 0; a < num; ++a) {
            size_t id = 0;
            double score = 0;
            iarc >> id;
            ids.push_back(id);
            iarc >> score;
            scores.push_back(score);
        }
    }
};

// Message Type
// store the vertex id with max score
struct score_message : public graphlab::IS_POD_TYPE {
    graphlab::vertex_id_type id;
    double score;
    score_message():id(0),score(0) { }
    score_message(graphlab::vertex_id_type id, double score) : id(id), score(score) { }

    score_message& operator += (const score_message& other) {
        id = other.id;
        score = other.score;
        return *this;
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << id << score;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> id >> score;
    }
};


// Find the children with max score for every valid vertex
class find_max_children : public graphlab::ivertex_program<graph_type, id_and_score, score_message>,
                          public graphlab::IS_POD_TYPE {
    // this is a local copy of the message                       
    double max_score;
    graphlab::vertex_id_type max_id;

public:
    void init(icontext_type& context, const vertex_type& vertex, const score_message& msg) { 
        max_score = msg.score;
        max_id = msg.id;
    }

    edge_dir_type gather_edges (icontext_type& context, const vertex_type& vertex) const {
        return graphlab:: OUT_EDGES;
    }

    id_and_score gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
         return id_and_score(edge.target().id(), edge.target().data().score);
    }

    //get values and make a vector
    void apply(icontext_type& context, vertex_type& vertex, const id_and_score& total) {
        std::cout<<"start to apply function"<<std::endl
        const std::vector<graphlab::vertex_id_type>& ids = total.ids;
        const std::vector<double>& scores = total.scores;
        max_score = scores[0];
        max_id = ids[0];
        if (ids.size()==0){
            std::cout << "The is Leaf: " <<std::endl;
            return;
        }
        for (size_t i = 1; i < ids.size(); ++i) {
            if (scores[i] < max_score) {
                max_score = scores[i];
                max_id = ids[i];
                std::cout << "MAX_ID: " << max_id << std::endl;
            }
        }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return graphlab::OUT_EDGES; 
    }

    void scatter (icontext_type& context, const vertex_type& vertex, edge_type& edge) const { 
        const score_message msg(max_score,max_id);
        if (edge.target().id() != msg.id) {
            edge.target().data().valid = false;
        }
        else{
            context.signal(edge.target(),msg);
        }        
    }
};


/** Return type for merge program
*/
struct return_vertex_merge_data
{
    std::vector<graphlab::vertex_id_type> ids;
    std::vector<double> scores;
    std::vector<>
    // std::vector<graphlab::vertex_type> vertices;
    // vertex(): vertices() { }
    vertex(graphlab::vertex_type vertex): vertices() {
        vertices.push_back(vertex);
        // scores.push_back(score);
    }
    vertex& operator += (const vertex& other) {
        for (size_t i = 0; i < other.vertices.size(); ++i) {
            vertices.push_back(other.vertices[i]);
            // scores.push_back(other.scores[i]);
        }
        return *this;
    }

    void save(graphlab::oarchive& oarc) const {
        size_t num = vertices.size();
        oarc << num;
        for (size_t a = 0; a < num; ++a) 
            oarc << vertices[a] << scores[a];
    }

    void load(graphlab::iarchive& iarc) {
        vertices.clear();
        vertices.clear();
        size_t num = 0;
        iarc >> num;
        for (size_t a = 0; a < num; ++a) {
            size_t vertex = 0;
            // double score = 0;
            iarc >> vertex;
            vertices.push_back(vertex);
            // iarc >> score;
            // scores.push_back(score);
        }
    }
};





/**
 * Merge the content from the next node in path to get final result
 */
class merge : public graphlab::ivertex_program<graph_type, id_and_score, score_message>,
                          public graphlab::IS_POD_TYPE {
    int offset;
    int offset_n;
    int length;
    int length_n;

public:
    // void init(icontext_type& context, const vertex_type& vertex, const score_message& msg) { 
    //     offset = msg.offset;
    //     offset_n = msg.offset_n;
    //     length = msg.length;
    //     length_n = msg.length_n;
    // }

    edge_dir_type gather_edges (icontext_type& context, const vertex_type& vertex) const {
        return graphlab:: OUT_EDGES;
    }


    id_and_score gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
         return vertex(edge.target().id(), edge.target().data().score);
    }

    //get values and make a vector
    void apply(icontext_type& context, vertex_type& vertex, const vertex& total) {
        // std::cout<<"start to apply function"<<std::endl
        const std::vector<graphlab::vertex>& vertices = total.vertices;
        // const std::vector<double>& scores = total.scores;
        // max_score = scores[0];
        // max_id = ids[0];
        offset = 0;
        offset_n = 0;
        length_n = 0;
        length = 0;
        substr = "";
        if (vertices.size()==0){
            std::cout << "The is Leaf: " <<std::endl;
            return;
        }
        for (size_t i = 1; i < vertices.size(); ++i) {
            if (vertices[i].data().vid == vertex.data().next_id) {
                offset = vertex.data().offset;
                length = vertex.data().length;
                offset_n = vertices[i].data().offset;
                length_n = vertices[i].data().offset;
                if (offset_n < (offset+lengh-1)){
                    vertex.data().merge_ct = vertex.data().content + 
                                           vertices[i].data().content.substr(offset+length-offset_n,length_n-offset-length+offset_n);    
                } 
            }
        }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return graphlab::OUT_EDGES; 
    }

    void scatter (icontext_type& context, const vertex_type& vertex, edge_type& edge) const { 
        const score_message msg(max_score,max_id);
        if (edge.target().id() != msg.id) {
            edge.target().data().valid = false;
        }
        else{
            context.signal(edge.target(),msg);
        }        
    }




        };





/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", save_vertex) to save the graph.
 */
struct best_path_writer {
    std::string save_vertex(const graph_type::vertex_type& vtx) {
        std::stringstream strm;
        strm << vtx.id() << "\t" << vtx.data().score << "\n";
        return strm.str();
    }
    std::string save_edge(graph_type::edge_type e) {
        std::stringstream strm;
        strm << "(" << e.source().id() << "\t" << e.target().id() << ")" << "\n";
        return strm.str();
    }
  //std::string save_edge(graph_type::edge_type e) { return ""; }
}; 


int main(int argc, char** argv) {

	// Initialize control plain using mpi
	graphlab::mpi_tools::init(argc, argv);
    graphlab::distributed_control dc;
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
    logger(LOG_INFO, "Assembly starting...\n");

    std::string infile;
    std::string outfile;
    graphlab::vertex_id_type source;
    graphlab::vertex_id_type destination;

    // Parse command line options -----------------------------------------------
    graphlab::command_line_options clopts("Welcome to assembly reads!");

    clopts.attach_option("infile", infile, "The aligned reads filename (required)");
    clopts.add_positional("infile");
    clopts.attach_option("outfile", outfile, "The filename for save the assembly results (required)");
    clopts.add_positional("outfile");

    // These two parameters should not be inputed by user in the future.
    clopts.attach_option("source", source, "The first reads of the sequence");
    clopts.attach_option("destination", destination, "The last reads of the sequence");
    clopts.print_description();

    std::cout << infile << " " << outfile << " " << source << " " << destination << std::endl;
    if(!clopts.parse(argc, argv)) {
    	dc.cout() << "Error in parsing command line arguments. \n";
    	return EXIT_FAILURE;
    }

    // Build the graph ----------------------------------------------------------
    graph_type graph(dc, clopts);

    if(!clopts.is_set("infile")) {
    	dc.cout() << "Input file not provided. \n";
    	return EXIT_FAILURE;
    }
    if(!clopts.is_set("outfile")) {
        dc.cout() << "Output file not provided. \n";
        return EXIT_FAILURE;
    }
    if(!clopts.is_set("source")) {
        dc.cout() << "The first read of the sequence not provided. \n";
        return EXIT_FAILURE;
    }
    if(!clopts.is_set("source")) {
        dc.cout() << "The last read of the sequence not provided. \n";
        return EXIT_FAILURE;
    }
    

    //graph.load(infile, line_parser);
    std::cout << infile << " " << outfile << " " << source << " " << destination << std::endl;
    if (!load_graph(infile, graph)) {
    	dc.cout() << "Cannot load the graph. \n";
    	return EXIT_FAILURE;
    }

    graph.finalize();

    // Running The Engine1 for exempt_reads_program -------------------------------------------------------
    graphlab::omni_engine<exempt_reads_program> engine1(dc, graph, "synchronous", clopts);
    engine1.signal(source, destination);
    engine1.start();
    const float runtime1 = engine1.elapsed_seconds();
    dc.cout() << "Finished Running engine1 in " << runtime1 << " seconds." << std::endl;

    // Running The Engine2 for find_max_children
        
    graphlab::omni_engine<find_max_children> engine2(dc, graph, "synchronous", clopts);
    engine2.signal(source, score_message(source, 0));
    engine2.start();
    const float runtime2 = engine2.elapsed_seconds();
    dc.cout() << "Finished Running engine2 in " << runtime2 << " seconds." << std::endl;

    // Save the final graph
    graph.save(outfile, best_path_writer(),
               false,    // do not gzip
               true,     // save vertices
               true);   // do not save edges

    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;

}

