
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <list>
#include <sstream>
#include <math.h>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <memory>
#include <cstring>

struct SiteBuffer {
  std::list<int> pos;
  std::list<int> n;
  std::list<int> k;
  std::list<int> n_mock;
  std::list<int> k_mock;
};

enum class ParsingState {
  SEQID, POSITION, NUCLEOTIDE, COVERAGE, COVERAGE_STRING, QUALITY_STRING, DONE
};

enum class CoverageState {
  NONE, FWD, REV, BEGIN, END
};

struct Site {
  int pos;
  int n;
  int k;
};

struct State {
  ParsingState parsing_state;
  CoverageState coverage_state;
  int remaining_coverage;
};

struct FullSite {
  int pos;
  int n_factor;
  int k_factor;
  int n_mock;
  int k_mock;

  bool operator==(const FullSite &other) const {
    return (n_factor == other.n_factor
            && k_factor == other.k_factor
            && n_mock == other.n_mock
            && k_mock == other.k_mock);
  }
};

namespace std {

  template <>
  struct hash<FullSite>
  {
    std::size_t operator()(const FullSite& k) const
    {
      size_t hash;
      hash = k.n_factor;
      hash *= 37;
      hash += k.k_factor;
      hash *= 37;
      hash += k.n_mock;
      hash *= 37;
      hash += k.k_mock;
      return hash;
    }
  };
}

constexpr int buffer_size = 100000;

bool parse_mpileup_line(std::array<char, buffer_size>& buffer, Site& plus_site, Site& minus_site, State& state);
inline void insert_site(SiteBuffer& buf, Site& factor_site, Site& mock_site);
inline void delete_site(SiteBuffer& buf);
inline void aggregate_weighted_sum(SiteBuffer& buf, std::vector<double>& weights, FullSite& output);
inline void count_site(std::unordered_map<FullSite, size_t>& map, FullSite& site);
inline void write_site(std::ofstream& file_stream, FullSite& site, std::string strand, int threshold);
inline void write_statistics(std::string& statistics_filename, std::unordered_map<FullSite, size_t>& statistics);
inline void reset_site(Site& site);
inline void reset_state(State& state);
inline int extract_integer(char* array, int start, int length);


struct Arguments {
  std::string region;
  std::string genome_fasta;
  std::string factor_bam;
  std::string mock_bam;
  std::string sites_file;
  std::string statistics_file;
  int window_size;
};


void gather_mockinbird_data(Arguments& args) {

  std::ofstream sites_file;
  sites_file.open(args.sites_file);

  std::vector<double> weights;
  int half_window = (args.window_size - 1) / 2;

  // at hoc decision - no math behind the heuristic
  double sigma2;
  if (half_window == 0) {
    sigma2 = 1;
  } else {
    sigma2 = half_window;
  }

  for (int i = 0; i < args.window_size; i++) {
    double weight = exp(- 0.5 * (i - half_window) * (i - half_window) / sigma2);
    weights.push_back(weight);
  }

  auto region = args.region;

  std::ostringstream factor_cmd_stream;
  factor_cmd_stream << "samtools mpileup -aa -C 0 -d 1000000000 -q 0 -Q 0 -f " << args.genome_fasta << " "
                    << args.factor_bam << " " << "-r '" << args.region << "'";
  auto cmd_factor = factor_cmd_stream.str();

  std::ostringstream mock_cmd_stream;
  mock_cmd_stream << "samtools mpileup -aa -C 0 -d 1000000000 -q 0 -Q 0 -f " << args.genome_fasta << " "
                    << args.mock_bam << " " << "-r '" << args.region << "'";
  auto cmd_mock = mock_cmd_stream.str();

  std::shared_ptr<FILE> factor_pipe(popen(cmd_factor.c_str(), "r"), pclose);
  if (!factor_pipe) throw std::runtime_error("factor read failed!");
  std::shared_ptr<FILE> mock_pipe(popen(cmd_mock.c_str(), "r"), pclose);
  if (!mock_pipe) throw std::runtime_error("mock read failed!");

  Site factor_site_plus;
  Site factor_site_minus;
  Site mock_site_plus;
  Site mock_site_minus;

  std::array<char, buffer_size> buffer;

  SiteBuffer plus_buffer;
  SiteBuffer minus_buffer;

  FullSite full_plus_site;
  FullSite full_minus_site;

  State factor_state;
  State mock_state;

  std::unordered_map<FullSite, size_t> statistics;

  unsigned n_sites = 0;
  while (!feof(factor_pipe.get()) & !feof(mock_pipe.get())) {

    reset_state(factor_state);
    reset_site(factor_site_plus);
    reset_site(factor_site_minus);
    if((fgets(buffer.data(), buffer_size, factor_pipe.get())) == nullptr) {
      continue;
    }
    while(parse_mpileup_line(buffer, factor_site_plus, factor_site_minus, factor_state)) {
      fgets(buffer.data(), buffer_size, factor_pipe.get());
    }

    reset_state(mock_state);
    reset_site(mock_site_plus);
    reset_site(mock_site_minus);
    if (fgets(buffer.data(), buffer_size, mock_pipe.get()) == nullptr) {
      continue;
    }
    while(parse_mpileup_line(buffer, mock_site_plus, mock_site_minus, mock_state)) {
      fgets(buffer.data(), buffer_size, mock_pipe.get());
    }

    insert_site(plus_buffer, factor_site_plus, mock_site_plus);
    insert_site(minus_buffer, factor_site_minus, mock_site_minus);
    n_sites += 1;

    // check if buffer is ready
    if (n_sites >= args.window_size) {
      aggregate_weighted_sum(plus_buffer, weights, full_plus_site);
      aggregate_weighted_sum(minus_buffer, weights, full_minus_site);

      count_site(statistics, full_plus_site);
      count_site(statistics, full_minus_site);

      write_site(sites_file, full_plus_site, "+", 1);
      write_site(sites_file, full_minus_site, "-", 1);

      delete_site(plus_buffer);
      delete_site(minus_buffer);
    }

  }
  write_statistics(args.statistics_file, statistics);
}

bool parse_mpileup_line(std::array<char, buffer_size>& buffer, Site& plus_site, Site& minus_site, State& state) {
  int i = 0;
  char ch;

  int pos_start = -1;
  int cov_start = -1;

  int line_done = false;

  while (!line_done && i < buffer_size) {
    ch = buffer[i];
    switch (state.parsing_state) {
      case ParsingState::SEQID:
        if (ch == '\t') {
          pos_start = i + 1;
          state.parsing_state = ParsingState::POSITION;
        };
        break;
      case ParsingState::POSITION:
        if (ch == '\t') {
          int pos_end = i - 1;
          int pos = extract_integer(buffer.data(), pos_start, pos_end);
          plus_site.pos = pos - 1;
          minus_site.pos = pos + 1;
          state.parsing_state = ParsingState::NUCLEOTIDE;
        };
        break;
      case ParsingState::NUCLEOTIDE:
        // we do not care
        cov_start = i + 2;
        i = cov_start;
        state.parsing_state = ParsingState::COVERAGE;
        break;
      case ParsingState::COVERAGE:
        if (ch == '\t') {
          state.parsing_state = ParsingState::COVERAGE_STRING;
          int cov_end = i - 1;
          state.remaining_coverage = extract_integer(buffer.data(), cov_start, cov_end);
        }
        break;
      case ParsingState::COVERAGE_STRING:
        if (ch == '\t') {
          state.parsing_state = ParsingState::QUALITY_STRING;
          break;
        }
        switch (ch) {
          case '^':
            i++; // always comes as ^~
            state.coverage_state = CoverageState::BEGIN;
            break;
          case 'A':
          case 'C':
          case 'G':
          case 'T':
          case '.':
            if (state.coverage_state == CoverageState::BEGIN) {
              ++plus_site.k;
            }
            ++plus_site.n;
            state.coverage_state = CoverageState::FWD;
            break;
          case '$':
            if (state.coverage_state == CoverageState::REV) {
              ++minus_site.k;
            }
            state.coverage_state = CoverageState::END;
            break;
          case 'a':
          case 'c':
          case 'g':
          case 't':
          case ',':
            ++minus_site.n;
            state.coverage_state = CoverageState::REV;
            break;
        }
        break;
      case ParsingState::QUALITY_STRING: {
        int remaining_slots = buffer_size - i - 1; // -1 because last slot in buffer is '\0'
        // we need one slot for '\n'
        if (state.remaining_coverage < remaining_slots) {
          state.parsing_state = ParsingState::DONE;
          line_done = true;
        } else {
          // all the remaining slots in the buffer are taken by the quality string
          state.remaining_coverage -= remaining_slots;
          line_done = true;
        }
        break;
      }
      case ParsingState::DONE:
        break;
    }
    i++;
  }

  return state.parsing_state != ParsingState::DONE;
}

inline void reset_site(Site& site) {
  site.n = 0;
  site.k = 0;
}

inline void reset_state(State& state) {
  state.parsing_state = ParsingState::SEQID;
  state.coverage_state = CoverageState::NONE;
}

void inline insert_site(SiteBuffer& buf, Site& factor_site, Site& mock_site) {
  buf.pos.push_back(factor_site.pos);
  buf.n.push_back(factor_site.n);
  buf.k.push_back(factor_site.k);
  buf.n_mock.push_back(mock_site.n);
  buf.k_mock.push_back(mock_site.k);
}

void inline delete_site(SiteBuffer& buf) {
  buf.pos.pop_front();
  buf.n.pop_front();
  buf.k.pop_front();
  buf.n_mock.pop_front();
  buf.k_mock.pop_front();
}

void inline count_site(std::unordered_map<FullSite, size_t>& map, FullSite& site) {
    if(map.find(site) != map.cend()) {
      ++map[site];
    } else {
      FullSite copy_site = site;
      ++map[copy_site];
    }
}

void inline write_site(std::ofstream& file_stream, FullSite& site, std::string strand, int threshold) {
  if (site.k_factor >= threshold) {
    file_stream << site.pos << '\t' << strand << '\t' << site.n_factor << '\t' << site.k_factor << '\t'
                << site.k_mock << std::endl;
  }
}

void inline aggregate_weighted_sum(SiteBuffer& buf, std::vector<double>& weights, FullSite& output) {
  auto window_size = weights.size();
  int mid_pos = window_size / 2;
  auto pos_it = buf.pos.cbegin();
  auto n_it = buf.n.cbegin();
  auto k_it = buf.k.cbegin();
  auto n_mock_it = buf.n_mock.cbegin();
  auto k_mock_it = buf.k_mock.cbegin();

  int pos;
  double n = 0;
  double k = 0;
  double n_mock = 0;
  double k_mock = 0;

  for (int i = 0; i < window_size; i++) {
    n = n + *n_it * weights[i];
    k = k + *k_it * weights[i];
    n_mock = n_mock + *n_mock_it * weights[i];
    k_mock = k_mock + *k_mock_it * weights[i];
    if (i == mid_pos) {
      pos = *pos_it;
    }
    n_it++;
    k_it++;
    n_mock_it++;
    k_mock_it++;
    pos_it++;

  }
  output.n_factor = round(n);
  output.k_factor = round(k);
  output.n_mock = round(n_mock);
  output.k_mock = round(k_mock);
  output.pos = pos;
}




void inline write_statistics(std::string& statistics_filename, std::unordered_map<FullSite, size_t>& statistics) {
  std::ofstream statistics_file;
  statistics_file.open(statistics_filename);
  for (std::pair<FullSite, size_t> element : statistics)
  {
    FullSite site = element.first;
    statistics_file << site.n_factor << '\t' << site.k_factor << '\t' << site.n_mock << '\t' << site.k_mock << '\t'
                    << element.second << std::endl;
  }
}

inline int extract_integer(char* array, int start, int end) {
  int length = end - start + 1;
  char substr[length + 1];
  substr[length] = '\0';
  strncpy(substr,  array + start, length);
  return atoi(substr);
}


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

  gather_mockinbird_data(args);
}