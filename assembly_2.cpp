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
    int next_id;
//    std::list<int> parent_list;
    bool valid;
    read_data(): length(0), content(""), offset(0), score(0), end(0), next_id(0), valid(false){ }
	explicit read_data(int length, std::string content, int offset, double score, int end, int next_id, bool valid) : length(length), content(content), offset(offset), score(score), end(end), valid(valid) { } 
    void save(graphlab::oarchive& oarc) const {
      oarc << length << content << offset << score << end << next_id << valid;
    }
    void load(graphlab::iarchive& iarc) {
      iarc >> length >> content >> offset >> score >> end >> next_id >> valid;
    }
};

// struct edge_data
// {
// 	double weight;
// 	edge_data(double weight) : weight(weight) { }
// };

// The type of graph used in this program
typedef graphlab::distributed_graph<read_data, graphlab::empty> graph_type;


//Create graph line by line
//Input file format: ReadID, Read Length, Read, Match Position, Match Score, <edges>
// bool line_parser(graph_type& graph, const std::string& filename, const std::string& textline) {

//     std::stringstream strm(textline);
//     graphlab::vertex_id_type vid;
//     int length;
//     std::string content;
//     int offset;
//     double score;

//     // first entry in the line is a vertex ID
//     strm >> vid;
//     strm >> length;
//     strm >> content;
//     strm >> offset;
//     strm >> score;
//     int end = offset + length;

//     // insert this read
//     graph.add_vertex(vid, read_data(length, content, offset, score, end, true));

//     // while there are elements in the line, continue to read until we fail
//     while(1){
//         graphlab::vertex_id_type other_vid;
//         strm >> other_vid;
//         if (strm.fail()) 
//             return false;
//         graph.add_edge(vid, other_vid);
//     }
//     return true;
// }

bool load_graph(const std::string& filename, graph_type& graph) {

	std::ifstream fin(filename.c_str());
    std::cout << "Load " << filename << " now ... " << std::endl;
	if(!fin.good()) return false;

	std::string textline;
	if (fin.is_open()) {
        std::cout << "File is open! " << std::endl;
		graphlab::vertex_id_type vid;
 //       graphlab::vertex_id_type last_vid = 0;
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

            // std::cout << "check last_vid..." << std::endl;
            // std::cout << "last_vid: " << last_vid << std::endl;

            // if (last_vid != 0){
            //     std::cout << "xxxxx!" <<std::endl;

            //     if(graph.num_vertices() > 0) {std::cout << last_vid << " not empty" << std::endl;}

            //     std::cout << graph.vertex(last_vid).id() << std::endl;
            //     std::cout << "last_vid is not zero! out edges = " << graph.vertex(last_vid).num_out_edges() << std::endl;
            //     if (graph.vertex(last_vid).num_out_edges() == 0)
            //         graph.add_edge(last_vid,vid);
            //     std::cout << "edge: (" << last_vid << ", " << vid << ")" << std::endl;
            // }
            // last_vid = vid;

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
        //std::cout<<"EXEMPT_INI"<<std::endl;
	last_id = msg;
    }

    edge_dir_type gather_edges (icontext_type& context, const vertex_type& vertex) const {
	//std::cout << "EXEMPT_GARTHER" <<std::endl;
        return graphlab:: NO_EDGES;
    }

    void apply(icontext_type& context, vertex_type& vertex, const graphlab::empty& empty){
        // if a vertex do not have a child and is not the last vertex, then tag it as an invalid vertex.
        //std::cout << "Exempt_apply" << std::endl;
        //std::cout << "Exempt_ID: " << vertex.id() << std::endl;
        if (vertex.num_out_edges() == 0 && vertex.id() != last_id) {
            vertex.data().valid = false;
        }
    }
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
	//std::cout << "EXEmpt_scatter"<< std::endl;
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
        //std::cout << "FIND_ININ" << std::endl;
        max_score = msg.score;
        max_id = msg.id;
    }

    edge_dir_type gather_edges (icontext_type& context, const vertex_type& vertex) const {
        //std::cout << "FIND_GARTHER" <<std::endl;
        return graphlab:: OUT_EDGES;
   }

    id_and_score gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
        //std::cout << "FIND_GARTHER2" <<std::endl;
        return id_and_score(edge.target().id(), edge.target().data().score);
    }

    //get values and make a vector
    void apply(icontext_type& context, vertex_type& vertex, const id_and_score& total) {

        std::cout << "The id now is: " << vertex.id() <<std::endl;
        const std::vector<graphlab::vertex_id_type>& ids = total.ids;
       // std::cout << "Here is wrong!" << std::endl;
        const std::vector<double>& scores = total.scores;
      //  std::cout << "Here 2" << std::endl;
       // std::cout << "socres.size is: " << scores.size() << std::endl;
        if(scores.size()==0){
           std::cout << "This is Leaf now " <<std::endl;
           return;
         }
        max_score = scores[0];
        //std::cout << "score[0] " << std::endl;
        max_id = ids[0];
        //std::cout << "FIND_Apply"<<std::endl;
   	//std::cout << "SIZE: " << ids.size() << std::endl;
       // if(ids.size()==0){
       //     std::cout << "This is Leaf now " <<std::endl;
       //    return;
       // }
        for (size_t i = 0; i < ids.size(); ++i) {
            std::cout << "find_here1" <<std::endl;
            std::cout << "score for " << i << "is: " << scores[i] << std::endl;
            std::cout << "max score is: " << max_score << std::endl; 
            if (scores[i] >  max_score) {
                std::cout << "find_HERE" <<std::endl;
                max_score = scores[i];
                max_id = ids[i];
                std::cout << "MAX_ID:" <<max_id <<std::endl;      
            }
        }
        vertex.data().next_id = max_id
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
        strm << vtx.id() << "\t" << vtx.data().score << vtx.data().next_id << "\n";
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
    dc.cout() << "Start to EXEMPT" <<std::endl;
    graphlab::omni_engine<exempt_reads_program> engine1(dc, graph, "synchronous", clopts);
    engine1.signal(source, destination);
    engine1.start();
    const float runtime1 = engine1.elapsed_seconds();
    dc.cout() << "Finished Running engine1 in " << runtime1 << " seconds." << std::endl;

    // Running The Engine2 for find_max_children
    dc.cout() << "Start to FIND" <<std::endl;
    graphlab::omni_engine<find_max_children> engine2(dc, graph, "synchronous", clopts);
    engine2.signal_all();  
 // engine2.signal(source, score_message(source, 0));
    engine2.start();
    const float runtime2 = engine2.elapsed_seconds();
    dc.cout() << "Finished Running engine2 in " << runtime2 << " seconds." << std::endl;

    // Save the final graph
    graph.save(outfile, best_path_writer(),
               false,    // do not gzip
               true,     // save vertices
          //     false);   // do not save edges
               true); //save edges
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;

}

