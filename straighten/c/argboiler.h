#include <argtable2.h>

#define ARG_CONSTS_(type,str,meth,lval,id,shrt,lng,desc,def) const type id##_default = def;
#define ARG_CONSTS(type,id,shrt,lng,desc,def) ARG_CONSTS_(type,id,shrt,lng,desc,def)
#define ARG_STRUCT_(type,str,meth,lval,id,shrt,lng,desc,def) type id;
#define ARG_STRUCT(type,id,shrt,lng,desc,def) ARG_STRUCT_(type,id,shrt,lng,desc,def)
#define ARG_TABVAL_(type,str,meth,lval,id,shrt,lng,desc,def) struct arg_##str * id = arg_##str##meth(shrt,lng,"<" #str ">", desc " (default " #def ")");
#define ARG_TABVAL(type,id,shrt,lng,desc,def) ARG_TABVAL_(type,id,shrt,lng,desc,def)
#define ARG_TABLE_(type,str,meth,lval,id,shrt,lng,desc,def) id,
#define ARG_TABLE(type,id,shrt,lng,desc,def) ARG_TABLE_(type,id,shrt,lng,desc,def)
#define ARG_TABINIT_(type,str,meth,lval,id,shrt,lng,desc,def) id->lval[0] = def;
#define ARG_TABINIT(type,id,shrt,lng,desc,def) ARG_TABINIT_(type,id,shrt,lng,desc,def)
#define ARG_INIT_(type,str,meth,lval,id,shrt,lng,desc,def) params->id = id->lval[0];
#define ARG_INIT(type,id,shrt,lng,desc,def) ARG_INIT_(type,id,shrt,lng,desc,def)

#define ARG_STR0 const char*,str,0,sval
#define ARG_STR1 const char*,str,1,sval
#define ARG_INT0 int,int,0,ival
#define ARG_INT1 int,int,1,ival
#define ARG_DBL0 double,dbl,0,dval
#define ARG_DBL1 double,dbl,1,dval

ARGBOILER(ARG_CONSTS)

typedef struct {
  ARGBOILER(ARG_STRUCT)
} args_t;

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
