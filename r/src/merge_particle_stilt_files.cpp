#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <iostream>
#include <assert.h>

// Heap item only stores sort keys and file index. The line payload is kept in
// per-file state to avoid copying large strings on every heap push/pop.
struct HeapItem {
    double timestamp;
    int id;
    std::size_t fileIndex;
};

struct HeapItemLess {
    bool operator()(const HeapItem& lhs, const HeapItem& rhs) const {
        // Priority queue places "largest" item at top. We want descending
        // timestamps, then ascending particle id.
        if (lhs.timestamp < rhs.timestamp) return true;
        if (lhs.timestamp > rhs.timestamp) return false;
        return lhs.id > rhs.id;
    }
};

struct FileState {
    std::ifstream stream;
    int offset = 0;
    double timestamp = 0.0;
    int id = 0;
    std::string rest;
};

static bool readNextRecord(FileState& file) {
    if (!(file.stream >> file.timestamp >> file.id)) {
        return false;
    }
    std::getline(file.stream, file.rest);
    file.id += file.offset;
    return true;
}


int main(int argc, char* argv[]) {
    // Check if there is at least one file to process
    if (argc < 4 || argc % 2 != 0) {
        std::cerr << "Usage: " << argv[0] << " <output_file_name> <file1> <offset1> <file2> <offset2> ..." << std::endl;
        return 1;
    }

    // Open all files and read first record from each.
    std::priority_queue<HeapItem, std::vector<HeapItem>, HeapItemLess> queue;
    std::vector<FileState> files;
    files.reserve((argc - 2) / 2);
    std::string output_filename = argv[1];

    std::string header;
    for (int i = 2; i < argc; i += 2) {
        files.emplace_back();
        FileState& file = files.back();
        file.stream.open(argv[i]);
        file.offset = std::stoi(argv[i + 1]);

        if (!file.stream.is_open()) {
            std::cerr << "Failed to open input file: " << argv[i] << std::endl;
            return 1;
        }

        // record the header
        if (!std::getline(file.stream, header)) {
            std::cerr << "Failed to read header from input file: " << argv[i] << std::endl;
            return 1;
        }
        //std::cout << "Read header: " << header << " from file " << argv[i] << std::endl;

        if (readNextRecord(file)) {
            queue.push({file.timestamp, file.id, files.size() - 1});
        }
    }

    // open output file
    std::ofstream out(output_filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << output_filename << std::endl;
        return 1;
    }

    // write the header recorded earlier
    out << header << '\n';

    while (!queue.empty()) {
        // get the highest-priority record (latest timestamp, then lowest id)
        HeapItem rec = queue.top();
        queue.pop();

        // write to output file
        FileState& file = files[rec.fileIndex];
        out << rec.timestamp << ' ' << rec.id << file.rest << '\n';

        // read the next record from the same file and add to the queue
        assert(rec.fileIndex < files.size());
        if (readNextRecord(file)) {
            queue.push({file.timestamp, file.id, rec.fileIndex});
        }
    }

    return 0;
}
