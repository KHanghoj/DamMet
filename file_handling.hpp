#pragma once

bool check_file_exists(std::string filename);

template <typename T>
void checkfilehandle(T &fh, std::string filename);

void filter_ref_sites(const std::string & selected_chrom, std::string & filename, char * ref);
void filter_ref_bed(const std::string & selected_chrom, std::string & filename, char * ref);


std::vector<std::string> parse_chrom_file(std::string & filename);
std::vector<std::pair<size_t, size_t>> parse_bed_file(std::string filename,  std::string &curr_chrom);
