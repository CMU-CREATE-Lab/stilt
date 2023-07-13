#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <assert.h>

// structure to hold information about each line
struct Record {
    double timestamp;
    int id;
    unsigned fileIndex;
    std::string rest;
    
    // Priority queue starts with highest value (e.g. descending order)

    bool operator<(const Record& other) const {
        // Sort primarily by timestamp by descending order since this is a backtrace.
        if (timestamp < other.timestamp) return true;
        if (timestamp > other.timestamp) return false;
        // Sort secondarily by particle id in ascending order (so we negate it)
        return -id < -other.id;
    }
};

int main(int argc, char* argv[]) {
    // Check if there is at least one file to process
    if (argc < 4 || argc % 2 != 0) {
        std::cerr << "Usage: " << argv[0] << " <output_file_name> <file1> <offset1> <file2> <offset2> ..." << std::endl;
        return 1;
    }

    // priority queue to hold the first record of each file
    std::priority_queue<Record> queue;

    // open all files and read the first record of each file
    std::vector<std::ifstream> files;
    std::vector<int> offsets;
    std::string output_filename = argv[1];

    std::string header;
    for (int i = 2; i < argc; i += 2) {
        files.emplace_back(argv[i]);
        offsets.push_back(std::stoi(argv[i + 1]));

        // record the header
        std::getline(files.back(), header);
        //std::cout << "Read header: " << header << " from file " << argv[i] << std::endl;

        // read the first record
        Record rec;
        files.back() >> rec.timestamp >> rec.id;
        std::getline(files.back(), rec.rest);
        rec.fileIndex = files.size() - 1;
        rec.id += offsets.back();
        queue.push(rec);
    }

    // open output file
    std::ofstream out(output_filename);

    // write the header recorded earlier
    out << header << '\n';

    while (!queue.empty()) {
        // get the record with smallest timestamp
        Record rec = queue.top();
        queue.pop();

        // write to output file
        out << rec.timestamp << ' ' << rec.id << rec.rest << '\n';

        // read the next record from the same file and add to the queue
        unsigned fileIndex = rec.fileIndex;
        assert(fileIndex < files.size());
        if (files[fileIndex] >> rec.timestamp >> rec.id) {
            std::getline(files[fileIndex], rec.rest);
            rec.id += offsets[fileIndex];
            queue.push(rec);
        }
    }

    return 0;
}
