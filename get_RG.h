// Parse the header, count the number of RG tags and return a list of their names
static bool count_RG(bam_hdr_t* hdr, size_t* count, char*** output_name)
{
  if (hdr->l_text < 3 ) {
    *count = 0;
    *output_name = NULL;
    return true;
  }
  kstring_t input = { 0, 0, NULL };
  kputsn(hdr->text, hdr->l_text, &input);

  //////////////////////////////////////////
  // First stage count number of @RG tags //
  //////////////////////////////////////////
  char* pointer = ks_str(&input);
  size_t n_rg = 0;
  // Guard against rare case where @RG is first header line
  // This shouldn't happen but could where @HD is omitted
  if (pointer[0] == '@' && pointer[1] == 'R' && pointer[2] == 'G' ) {
    ++n_rg;
    pointer += 3;
  }
  char* line;
  while ((line = strstr(pointer, "\n@RG")) != NULL) {
    ++n_rg;
    pointer = line + 1;
  }

  //////////////////////////////////
  // Second stage locate @RG ID's //
  //////////////////////////////////
  char** names = (char**)calloc(sizeof(char*), n_rg);
  size_t next = 0;

  regex_t rg_finder;
  if (regcomp(&rg_finder, "^@RG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE) != 0) {
    free(input.s);
    free(names);
    return false;
  }
  regmatch_t* matches = (regmatch_t*)calloc(sizeof(regmatch_t),2);
  int error;
  char* begin = ks_str(&input);

  while ((error = regexec(&rg_finder, begin, 2, matches, 0)) == 0) {
    kstring_t str = { 0, 0, NULL };
    kputsn(begin+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &str);
    names[next++] = ks_release(&str);
    begin += matches[0].rm_eo;
  }

  if (error != REG_NOMATCH) {
    // cleanup
    regfree(&rg_finder);
    free(matches);
    free(names);
    free(input.s);
    return false;
  }
  free(matches);

  // return results
  *count = n_rg;
  *output_name = names;
  regfree(&rg_finder);
  free(input.s);
  return true;
}
