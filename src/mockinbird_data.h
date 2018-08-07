//
// Created by Christian Roth on 7/31/18.
//

#ifndef MOCKINBIRD_MOCKINBIRD_DATA_H
#define MOCKINBIRD_MOCKINBIRD_DATA_H

#include <list>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <iostream>
#include <fstream>

constexpr size_t ASCII_len = 128;
constexpr int buffer_size = 2000000;

struct SiteBuffer {
  std::list<int> pos;
  std::list<int> n;
  std::list<int> k;
  std::list<int> n_mock;
  std::list<int> k_mock;
};

enum class ParsingState {
  SEQID, POSITION, NUCLEOTIDE, COVERAGE, COVERAGE_STRING
};

struct State {
  ParsingState parsing_state;
  bool skip;
  std::array<unsigned int, ASCII_len> coverage_counts;
};

inline void reset_state(State& state) {
  state.parsing_state = ParsingState::SEQID;
  state.skip = false;
  state.coverage_counts['A'] = 0;
  state.coverage_counts['C'] = 0;
  state.coverage_counts['G'] = 0;
  state.coverage_counts['T'] = 0;
  state.coverage_counts['.'] = 0;
  state.coverage_counts['a'] = 0;
  state.coverage_counts['c'] = 0;
  state.coverage_counts['g'] = 0;
  state.coverage_counts['t'] = 0;
  state.coverage_counts[','] = 0;
}

struct Site {
  int pos;
  bool plus_strand;
  int n;
  int k;
};




enum class CoverageState {
  NONE, FWD, REV, BEGIN, END
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

inline char reverse_complement(char nucleotide) {
  char compl_nuc;
  switch(nucleotide) {
    case 'A':
      compl_nuc = 'T';
      break;
    case 'C':
      compl_nuc = 'G';
      break;
    case 'G':
      compl_nuc = 'C';
      break;
    case 'T':
      compl_nuc = 'A';
      break;
  }
  return compl_nuc;
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

void inline write_site(std::ofstream& file_stream, FullSite& site, char strand, int threshold) {
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

struct ParseData {
  std::array<char, buffer_size> buffer;
  char ref_nucleotide_plus;
  char mut_nucleotide_plus;
  char ref_nucleotide_minus;
  char mut_nucleotide_minus;
  State parsing_state;
};


void parse_mpileup_line(ParseData& data, Site& site) {

  auto& ref_nuc_plus = data.ref_nucleotide_plus;
  auto& mut_nuc_plus = data.mut_nucleotide_plus;
  auto& ref_nuc_minus = data.ref_nucleotide_minus;
  auto& mut_nuc_minus = data.mut_nucleotide_minus;
  auto& state =  data.parsing_state;
  auto& buffer = data.buffer;

  int i = 0;
  char ch;

  int pos_start = -1;

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
          site.pos = pos;
          state.parsing_state = ParsingState::NUCLEOTIDE;
        };
        break;
      case ParsingState::NUCLEOTIDE:
        if (ch == ref_nuc_plus) {
          site.plus_strand = true;
        } else if (ch == ref_nuc_minus) {
          site.plus_strand = false;
        } else {
          state.skip = true;
          line_done = true;
        }
        ++i;
        state.parsing_state = ParsingState::COVERAGE;
        break;
      case ParsingState::COVERAGE:
        if (ch == '\t') {
          state.parsing_state = ParsingState::COVERAGE_STRING;
        }
        break;
      case ParsingState::COVERAGE_STRING:
        if (ch == '\t') {
          line_done = true;
          break;
        } else {
          state.coverage_counts[ch] += 1;
        }
        break;
    }
    i++;
  }

  if (site.plus_strand) {
    site.n = (
      state.coverage_counts['.']
      + state.coverage_counts['A']
      + state.coverage_counts['C']
      + state.coverage_counts['G']
      + state.coverage_counts['T']
    );
    site.k = state.coverage_counts['C'];
  } else {
    site.n = (
      state.coverage_counts[',']
      + state.coverage_counts['a']
      + state.coverage_counts['c']
      + state.coverage_counts['g']
      + state.coverage_counts['t']
    );
    site.k = state.coverage_counts['g'];
  }
}

#endif //MOCKINBIRD_MOCKINBIRD_DATA_H
