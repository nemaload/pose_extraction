#include "argboiler_args.h"
#include "argboiler.h"
#include <stdlib.h>

struct arg_lit* arg_lit_n(const char* shortopts, const char* longopts, const char* datatype, const char* glossary) {
  return arg_litn(shortopts,longopts,0,10,glossary);
}

struct arg_lit* arg_lit_0(const char* shortopts, const char* longopts, const char* datatype, const char* glossary) {
  return arg_lit0(shortopts,longopts,glossary);
}

ARGBOILER(ARG_CONSTS)

int parse_args(int argc, char** argv, args_t* params) {
  /*
   * We use the argtable library to parse arguments. The first step is
   * building the titular argtable. This is the third place.
   */
  ARGBOILER(ARG_TABVAL)
  struct arg_end* end = arg_end(10);
  /*
   * Don't forget to actually add your new args into the argtable below
   * before the "end" sentinel.
   */
  void* argtable[] = {ARGBOILER(ARG_TABLE) end};
  int nerrors;

  /*
   * Here's the fifth place, where we fill in the defaults prior to
   * running the arg_parse function.
   */

  ARGBOILER(ARG_TABINIT)

  nerrors = arg_parse(argc, argv, argtable);
  if(nerrors > 0) {
    printf("\e[A");
    arg_print_errors(stderr, end, argv[0]);
    fprintf(stderr,"\nUsage:\n%s",argv[0]);
    arg_print_syntaxv(stderr, argtable, "\n\n");
    arg_print_glossary(stderr, argtable, "\t%-25s %s\n");
    fprintf(stderr,"\n");
    exit(nerrors);
  } else {
    /*
     * Finally, here we fill in the "params" structure.
     *
     * TODO:
     * There might be some way to handle these five needs
     * while specifying arguments in only one place, by
     * some very clever macro trickery, involving using a
     * functional macro name as the argument to a functional
     * macro, and defining such a higher-order macro with
     * all the param data. I'll get to that later...
     */
    ARGBOILER(ARG_INIT)
    return 0;
  }
}
