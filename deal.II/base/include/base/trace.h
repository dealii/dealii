/* Tracing macros for logging purpose if DEBUG is set. */
#ifndef __trace_h
#define __trace_h
#include <base/logstream.h>
#ifdef DEBUG
// Trace a variable, i.e. print variable name and value to deallog
#define TRACEVAR(arg) deallog << __LINE__ << ": " << # arg << " = " << (arg) << endl; 
// Trace a message, i.e. print message to deallog
#define TRACEMSG(arg) deallog << __LINE__ << ": " << arg << endl;
// Trace vector, i.e. print vector components to deallog
#define TRACEVEC(arg) deallog << __LINE__ << ": " << #arg << " = "; \
                      for(int TRACE_i = 0; TRACE_i < arg.size(); TRACE_i ++) deallog << arg(TRACE_i) << " "; \
                      deallog << endl;
#else
#define TRACEVAR(arg)
#define TRACEMSG(arg)
#define TRACEVEC(arg)
#endif
#endif
