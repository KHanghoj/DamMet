#pragma once

faidx_t * ref_init(std::string & filename);
char * fetch_chrom(const faidx_t *fai, std::string & chrom);
