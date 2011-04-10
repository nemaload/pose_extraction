#ifndef ARGBOILER_
#include <argtable2.h>
#include <stdlib.h>

#define ARG_CONSTS_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) type id##_default;
#define ARG_CONSTS(type,id,shrt,lng,desc,def) ARG_CONSTS_(type,id,shrt,lng,desc,def,type##_DEF(def))
#define ARG_STRUCT_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) type id;
#define ARG_STRUCT(type,id,shrt,lng,desc,def) ARG_STRUCT_(type,id,shrt,lng,desc,def,type##_DEF(def))
#define ARG_TABVAL_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) struct arg_##strct * id = arg_##strct##meth(shrt,lng,typestr, desc defstr);
#define ARG_TABVAL(type,id,shrt,lng,desc,def) ARG_TABVAL_(type,id,shrt,lng,desc,def,type##_DEF(def))
#define ARG_TABLE_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) id,
#define ARG_TABLE(type,id,shrt,lng,desc,def) ARG_TABLE_(type,id,shrt,lng,desc,def,type##_DEF(def))
#define ARG_TABINIT_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) id->lval = id##_default = def;
#define ARG_TABINIT(type,id,shrt,lng,desc,def) ARG_TABINIT_(type,id,shrt,lng,desc,def,type##_DEF(def))
#define ARG_INIT_(type,strct,typestr,meth,lval,id,shrt,lng,desc,def,defstr) params->id = id->lval;
#define ARG_INIT(type,id,shrt,lng,desc,def) ARG_INIT_(type,id,shrt,lng,desc,def,type##_DEF(def))

#define ARG_FIL0 const char*,str,"<file>",0,sval[0]
#define ARG_FIL0_DEF(d) ""
#define ARG_FIL1 const char*,str,"<file>",1,sval[0]
#define ARG_FIL1_DEF(d) ""
#define ARG_STR0 const char*,str,"<string>",0,sval[0]
#define ARG_STR0_DEF(d) " [default: " d "]"
#define ARG_STR1 const char*,str,"<string>",1,sval[0]
#define ARG_STR1_DEF(d) ""
#define ARG_INT0 int,int,"<int>",0,ival[0]
#define ARG_INT0_DEF(d) " [default: " #d "]"
#define ARG_INT1 int,int,"<int>",1,ival[0]
#define ARG_INT1_DEF(d) ""
#define ARG_DBL0 double,dbl,"<real>",0,dval[0]
#define ARG_DBL0_DEF(d) " [default: " #d "]"
#define ARG_DBL1 double,dbl,"<real>",1,dval[0]
#define ARG_DBL1_DEF(d) ""
#define ARG_LIT0 int,lit,"",_0,count
#define ARG_LIT0_DEF(d) ""
#define ARG_LITN int,lit,"",_n,count
#define ARG_LITN_DEF(d) ""

#define ARGBOILER_

typedef struct {
  ARGBOILER(ARG_STRUCT)
} args_t;

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
   * building the titular argtable.
   */
  ARGBOILER(ARG_TABVAL)
  struct arg_end* end = arg_end(10);
  /*
   * Initialize argtable
   */
  void* argtable[] = {ARGBOILER(ARG_TABLE) end};
  int nerrors;

  /*
   * Fill in defaults
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
     */
    ARGBOILER(ARG_INIT)
    return 0;
  }
}

#endif
