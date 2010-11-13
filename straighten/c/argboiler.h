#ifndef ARGBOILER_
#include <argtable2.h>

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

int parse_args(int argc, char** argv, args_t* params);

#endif
