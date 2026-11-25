#include "output.h"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <stdexcept>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;

// add the deconvolution hyphen to the end of the BC or BX, wherever it is in the read
static std::string extend_barcode(std::string header, int extension)
{
    const char* whitespace = " \t\r\n\v\f";
    auto field = header.find_first_of(whitespace);

    while (field != std::string::npos) {
        field = header.find_first_not_of(whitespace, field);
        if (field == std::string::npos) {
            break;
        }

        auto end = header.find_first_of(whitespace, field);
        if (end == std::string::npos) {
            end = header.size();
        }

        if (end - field > 5 &&
            (header.compare(field, 5, "BX:Z:") == 0 ||
             header.compare(field, 5, "BC:Z:") == 0)) {
            header.insert(end, "-" + std::to_string(extension));
            break;
        }

        field = end;
    }

    return header;
}

void output(
    std::string inputFile,
    std::string outputFile,
    std::vector<Read>& reads,
    int min_length)
{
    const char format = sequence_format(inputFile);

    std::ifstream input_file;
    boost::iostreams::filtering_istream in;
    open_sequence_input(inputFile, input_file, in);

    std::ofstream output_file;
    boost::iostreams::filtering_ostream out;
    open_sequence_output(outputFile, output_file, out);

    std::size_t next_read_id = 0;
    std::string previous_key;
    std::string tag;
    std::string line;
    std::vector<std::string> buffer;

    const auto write_record = [&]() {
        if (!to_deconvolve(buffer, format, min_length, tag)) {
            return;
        }

        const auto key = mate_key(buffer[0], format);
        std::size_t read_id;
        if (next_read_id > 0 && key == previous_key) {
            read_id = next_read_id - 1;
        } else {
            read_id = next_read_id++;
            previous_key = key;
        }

        const auto extension = reads[read_id].barcode_extension;
        out << extend_barcode(buffer[0], extension)
            << '\t' << extension << '\n';
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == format) {
            write_record();
            buffer = {line};
        } else {
            buffer.push_back(line);
        }
    }

    write_record();
    finish_sequence_output(outputFile, output_file, out);
}

void fasta_output(
    std::string inputFile,
    std::string outputFile,
    std::vector<Read>& reads,
    int min_length)
{
    const char format = sequence_format(inputFile);

    std::ifstream input_file;
    boost::iostreams::filtering_istream in;
    open_sequence_input(inputFile, input_file, in);

    std::ofstream output_file;
    boost::iostreams::filtering_ostream out;
    open_sequence_output(outputFile, output_file, out);

    std::size_t next_read_id = 0;
    std::string previous_key;
    std::string tag;
    std::string line;
    std::vector<std::string> buffer;

    const auto write_record = [&]() {
        if (buffer.empty()) {
            return;
        }

        if (to_deconvolve(buffer, format, min_length, tag)) {
            const auto key = mate_key(buffer[0], format);
            std::size_t read_id;
            if (next_read_id > 0 && key == previous_key) {
                read_id = next_read_id - 1;
            } else {
                read_id = next_read_id++;
                previous_key = key;
            }

            const auto extension = reads[read_id].barcode_extension;
            out << extend_barcode(buffer[0], extension) << '\n';
        } else {
            out << buffer[0] << '\n';
        }

        for (std::size_t i = 1; i < buffer.size(); ++i) {
            out << buffer[i] << '\n';
        }
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == format) {
            write_record();
            buffer = {line};
        } else {
            buffer.push_back(line);
        }
    }

    write_record();
    finish_sequence_output(outputFile, output_file, out);
}
