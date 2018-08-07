
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
  std::string mock_bam;
  std::string sites_file;
  std::string statistics_file;
  char ref_nucleotide;
  char mut_nucleotide;
  int window_size;
};

class AggBuffer {

  unsigned size;
  unsigned buffer_head;

  std::vector<unsigned> pos;
  std::vector<unsigned> n;
  std::vector<unsigned> k;
  std::vector<unsigned> n_mock;
  std::vector<unsigned> k_mock;

  unsigned n_agg;
  unsigned k_agg;
  unsigned n_mock_agg;
  unsigned k_mock_agg;

public:
  AggBuffer(int size);
  void insert(Site& factor_site, Site& mock_site);
  void aggregate_site(FullSite& site);
  void discard_site();
};

AggBuffer::AggBuffer(int size) {
  this->size = size;
  pos.resize(size);
  n.resize(size);
  k.resize(size);
  n_mock.resize(size);
  k_mock.resize(size);
  buffer_head = 0;
}

void AggBuffer::insert(Site &factor_site, Site &mock_site) {
  int insert_pos = buffer_head - 1 + size;
  insert_pos %= size;

  n[insert_pos] = factor_site.n;
  k[insert_pos] = factor_site.k;
  n_mock[insert_pos] = mock_site.n;
  k_mock[insert_pos] = mock_site.k;
  pos[insert_pos] = factor_site.pos;

  n_agg += factor_site.n;
  k_agg += factor_site.k;
  n_mock_agg += mock_site.n;
  k_mock_agg += mock_site.k;
}

void AggBuffer::aggregate_site(FullSite &site) {
  unsigned mid_pos = (buffer_head + (size / 2)) % size;
  site.pos = pos[mid_pos];
  site.n_factor = n_agg;
  site.k_factor = k_agg;
  site.n_mock = n_mock_agg;
  site.k_mock = k_mock_agg;
}

void AggBuffer::discard_site() {
  n_agg -= n[buffer_head];
  k_agg -= k[buffer_head];
  n_mock_agg -= n_mock[buffer_head];
  k_mock_agg -= k_mock[buffer_head];
  buffer_head += 1;
  buffer_head %= size;

}

void gather_mockinbird_data(Arguments& args);
void parse_mpileup_line(ParseData& data, Site& site);


int main(int argc, char *argv[]) {

  Arguments args;
  int i = 0;
  args.factor_bam = std::string(argv[++i]);
  args.mock_bam = std::string(argv[++i]);
  args.genome_fasta = std::string(argv[++i]);
  args.region = std::string(argv[++i]);
  args.window_size = atoi(argv[++i]);
  args.sites_file = std::string(argv[++i]);
  args.statistics_file = std::string(argv[++i]);
  args.ref_nucleotide = *argv[++i];
  args.mut_nucleotide = *argv[++i];

  gather_mockinbird_data(args);
}

void gather_mockinbird_data(Arguments& args) {

  std::ofstream sites_file;
  sites_file.open(args.sites_file);

  auto region = args.region;

  std::ostringstream factor_cmd_stream;
  factor_cmd_stream << "samtools mpileup -aa -C 0 -d 500000 -q 0 -Q 0 -f " << args.genome_fasta << " "
                    << args.factor_bam << " " << "-r '" << args.region << "'";
  auto cmd_factor = factor_cmd_stream.str();

  std::ostringstream mock_cmd_stream;
  mock_cmd_stream << "samtools mpileup -aa -C 0 -d 500000 -q 0 -Q 0 -f " << args.genome_fasta << " "
                    << args.mock_bam << " " << "-r '" << args.region << "'";
  auto cmd_mock = mock_cmd_stream.str();

  std::shared_ptr<FILE> factor_pipe(popen(cmd_factor.c_str(), "r"), pclose);
  if (!factor_pipe) throw std::runtime_error("factor read failed!");
  std::shared_ptr<FILE> mock_pipe(popen(cmd_mock.c_str(), "r"), pclose);
  if (!mock_pipe) throw std::runtime_error("mock read failed!");

  ParseData data = ParseData{};
  data.ref_nucleotide_plus = args.ref_nucleotide;
  data.ref_nucleotide_minus = reverse_complement(args.ref_nucleotide);
  data.mut_nucleotide_plus = args.mut_nucleotide;
  data.mut_nucleotide_minus = tolower(args.mut_nucleotide);

  Site mock_site{};
  Site factor_site{};
  AggBuffer plus_buffer = AggBuffer{args.window_size};
  AggBuffer minus_buffer = AggBuffer{args.window_size};

  FullSite full_site;
  std::unordered_map<FullSite, size_t> statistics;

  unsigned n_sites = 0;
  while (!feof(factor_pipe.get()) & !feof(mock_pipe.get())) {

    reset_state(data.parsing_state);
    if ((fgets(data.buffer.data(), buffer_size, factor_pipe.get())) == nullptr) {
      continue;
    }
    parse_mpileup_line(data, factor_site);

    if (fgets(data.buffer.data(), buffer_size, mock_pipe.get()) == nullptr) {
      continue;
    }

    if (!data.parsing_state.skip) {
      reset_state(data.parsing_state);
      parse_mpileup_line(data, mock_site);

      AggBuffer& buffer = minus_buffer;
      char strand_symbol = '-';
      if (factor_site.plus_strand) {
        buffer = plus_buffer;
        strand_symbol = '+';
      }
      buffer.insert(factor_site, mock_site);

      n_sites += 1;

      // check if buffer is ready
      if (n_sites >= args.window_size) {
        buffer.aggregate_site(full_site);

        count_site(statistics, full_site);
        write_site(sites_file, full_site, strand_symbol, 1);

        buffer.discard_site();
      }
    }

  }
  write_statistics(args.statistics_file, statistics);
}