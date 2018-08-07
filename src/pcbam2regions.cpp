
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <list>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cstring>
#include "mockinbird_data.h"


struct Arguments {
  std::string region;
  std::string genome_fasta;
  std::string factor_bam;
  std::string region_file;
  char ref_nucleotide;
  char mut_nucleotide;
  int window_size;
};

struct Region {
  bool is_plus_strand;
  unsigned start;
  unsigned end;
  unsigned n_interesting_sites;
  bool is_interesting;
};

inline void write_region(std::ofstream& file_stream, Region& region) {
  char strand = '-';
  if (region.is_plus_strand) {
    strand = '+';
  }
  file_stream << region.start - 1 << '\t' << region.end << '\t' << strand << '\t' << std::endl;
}

class AggBuffer {

  unsigned size;
  unsigned buffer_head;
  unsigned buffer_tail;
  std::vector<unsigned> pos;
  std::vector<bool> interesting_sites;
  unsigned n_interesting_sites;

public:
  AggBuffer(int size);
  void insert(Site& factor_site);
  void aggregate_site(Region& region);
  void discard_site();
  unsigned n_processed_sites;
};

AggBuffer::AggBuffer(int size) {
  this->size = size;
  pos.resize(size);
  interesting_sites.resize(size);
  n_interesting_sites = 0;
  buffer_head = size - 1;
  buffer_tail = size - 1;
  n_processed_sites = 0;
}

void AggBuffer::insert(Site &site) {
  int insert_pos = (buffer_head + 1) % size;
  interesting_sites[insert_pos] = site.k > 0;
  pos[insert_pos] = site.pos;
  n_interesting_sites += interesting_sites[insert_pos];
  buffer_head = insert_pos;
  n_processed_sites += 1;

  if (n_processed_sites >= size) {
    buffer_tail = (buffer_tail + 1) % size;
  }
}

void AggBuffer::aggregate_site(Region &region) {
  unsigned mid_pos = (buffer_tail + (size / 2)) % size;
  region.is_interesting = interesting_sites[mid_pos];
  region.start = pos[buffer_tail];
  region.end = pos[buffer_head];
  region.n_interesting_sites = n_interesting_sites;
}

void AggBuffer::discard_site() {
  n_interesting_sites -= interesting_sites[buffer_tail];
}

void gather_mockinbird_data(Arguments& args);

int main(int argc, char *argv[]) {

  Arguments args;
  int i = 0;
  args.factor_bam = std::string(argv[++i]);
  args.genome_fasta = std::string(argv[++i]);
  args.region = std::string(argv[++i]);
  args.window_size = atoi(argv[++i]);
  args.region_file = std::string(argv[++i]);
  args.ref_nucleotide = *argv[++i];
  args.mut_nucleotide = *argv[++i];

  gather_mockinbird_data(args);
}

void gather_mockinbird_data(Arguments& args) {

  std::ofstream region_file;
  region_file.open(args.region_file);

  std::ostringstream factor_cmd_stream;
  factor_cmd_stream << "samtools mpileup -aa -C 0 -d 500000 -q 0 -Q 0 -f " << args.genome_fasta << " "
                    << args.factor_bam << " " << "-r '" << args.region << "'";
  auto cmd_factor = factor_cmd_stream.str();

  std::shared_ptr<FILE> factor_pipe(popen(cmd_factor.c_str(), "r"), pclose);
  if (!factor_pipe) throw std::runtime_error("factor read failed!");

  ParseData data = ParseData{};
  data.ref_nucleotide_plus = args.ref_nucleotide;
  data.ref_nucleotide_minus = reverse_complement(args.ref_nucleotide);
  data.mut_nucleotide_plus = args.mut_nucleotide;
  data.mut_nucleotide_minus = tolower(args.mut_nucleotide);

  Site site{};
  AggBuffer plus_buffer = AggBuffer{args.window_size};
  AggBuffer minus_buffer = AggBuffer{args.window_size};

  Region region;

  while (!feof(factor_pipe.get())) {

    reset_state(data.parsing_state);
    if ((fgets(data.buffer.data(), buffer_size, factor_pipe.get())) == nullptr) {
      continue;
    }
    parse_mpileup_line(data, site);
    region.is_plus_strand = site.plus_strand;

    if (!data.parsing_state.skip) {

      AggBuffer* buffer = &minus_buffer;
      if (region.is_plus_strand) {
        buffer = &plus_buffer;
      }
      buffer->insert(site);

      // check if buffer is ready
      if (buffer->n_processed_sites >= args.window_size) {

        buffer->aggregate_site(region);
        if (region.is_interesting && region.n_interesting_sites > 1) {
          write_region(region_file, region);
        }
        buffer->discard_site();
      }
    }
  }
}

