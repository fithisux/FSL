// 
// Implementation of a utility that will allow you
// to get information similar to bt on gdb if/when 
// your program crashes, even when not compiling with
// debug information. All you need to do is
//
// 1. Put 
//    #include "stack_dump.h"
//    among your includes.
// 2. Include the line
//    StackDump::Install()
//    at the start of main()
// 3. Add -rdynamic to your USRLDFLAGS
// 4. Wait for your program to crash
//
// Caveats:
// 1. When compiling with all optimisations some
//    functions will be inlined and will not appear
//    on the stack.
// 2. Highly non-portable. It works on the current systems 
//    on my Macbook, jalapeno and jalapeno00. Anything else 
//    is a bonus.
// 
// stack_dump.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

#ifndef stack_dump_h
#define stack_dump_h

#include <iostream>
#include <string>
#include <cstdio>
#include <csignal>
#include <execinfo.h>
#include <cxxabi.h>

class StackDump
{
private:
  static const unsigned int MAX_SDEPTH = 256;
  static void print_backtrace(int sig, siginfo_t *info, void *dodgy);
public:
  static void Install();
};

void StackDump::Install()
{
  struct sigaction sa;
  sa.sa_sigaction = print_backtrace;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO;
  sigaction(SIGSEGV,&sa,NULL);
  sigaction(SIGBUS,&sa,NULL);
  sigaction(SIGFPE,&sa,NULL);
}

#ifdef __APPLE__
void StackDump::print_backtrace(int sig, siginfo_t *info, void *dodgy)
{
  // Start by reporting on nature of signal
  std::string nature;
  if (sig == SIGSEGV) {
    nature = "Segmentation violation";
    switch (info->si_code) {
    case SEGV_MAPERR:
      nature += ", Address not mapped"; break;
    case SEGV_ACCERR:
      nature += ", Invalid permission for address"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Offending address = %p",info->si_addr);
    nature += tmp;
  }
  else if (sig == SIGBUS) {
    nature = "Bus error";
    switch (info->si_code) {
    case BUS_ADRALN:
      nature += ", Invalid address alignment"; break;
    case BUS_ADRERR:
      nature += ", Non-existent physical address"; break;
    case BUS_OBJERR:
      nature += ", Object specific hardware error"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Offending address = %p",info->si_addr);
    nature += tmp;
  }
  else if (sig == SIGFPE) {
    switch (info->si_code) {
    case FPE_INTDIV:
      nature += ", Integer divide-by-zero"; break;
    case FPE_INTOVF:
      nature += ", Integer overflow"; break;
    case FPE_FLTDIV:
      nature += ", Floating point divide-by-zero"; break;
    case FPE_FLTOVF:
      nature += ", Floating point overflow"; break;
    case FPE_FLTUND:
      nature += ", Floating point underflow"; break;
    case FPE_FLTRES:
      nature += ", Floating point loss of precision"; break;
    case FPE_FLTINV:
      nature += ", Invalid floating point operation"; break;
    case FPE_FLTSUB:
      nature += ", Subscript out of range"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Address of offending operation = %p",info->si_addr);
    nature += tmp;
  }
  std::cout << nature << std::endl;

  // Unwind stack
  void *stack_addr[MAX_SDEPTH];
  char **stack_strings;
  size_t stack_depth = backtrace(stack_addr,MAX_SDEPTH);
  // It seems on MacOS we have lost the address to the
  // function that caused the signal. The following, quite
  // dodgy, code fixes that. From my tests on jalapeno and
  // jalapeno00 the address is still there so we don't need
  // this hack for Linux (fingers crossed).
  ucontext_t *uc = (ucontext_t *) dodgy;
#ifdef __i386__  // If compiling for a 32-bit cpu
  stack_addr[2] = (void *) uc->uc_mcontext->__ss.__eip;
#else // If compiling for 64-bit cpu, different register names compared to 32-bit 
  stack_addr[2] = (void *) uc->uc_mcontext->__ss.__rip;
#endif 
  stack_strings = backtrace_symbols(stack_addr,stack_depth);

  // Demangle C++ function names. Prior to that the strings from
  // backtrace_symbols have to be "cleaned up", a non-portable
  // and empirical abomination.  
  // The "cleaning up" of the strings below is based on empirical 
  // tests on my 64-bit i7 Macbook running Snowleopard.   
  for (unsigned int i=2; i<stack_depth; i++) {
    char bls[256]; strcpy(bls,stack_strings[i]);
    char *file = NULL; char *fn=NULL; char *fend=NULL;
    for (file=bls; *file!=' '; file++)      //
      ;                                     //
    for (; *file==' '; file++)              // Skip confusing frame #
      ;                                     //
    for (fn=file; !(*fn=='0' && *(fn+1)=='x'); fn++)
      ; // We have found address
    for (; *fn!=' '; fn++)
      ; *fn='\0'; fn++; // Now at start of function name
    for (fend=fn; *fend!=' '; fend++)
      ; *fend='\0';
    size_t length; int status;
    char *demangled = abi::__cxa_demangle(fn,NULL,&length,&status);
    if (demangled) {
      std::cout << std::string(file) << " " << std::string(demangled) << std::endl;
      free(demangled);
    }
    else std::cout << std::string(file) << " " << std::string(fn) << std::endl; 
  }
  
  free(stack_strings);

  exit(EXIT_FAILURE);
  return;
}
#endif // #ifdef __APPLE__

#ifdef __linux__
void StackDump::print_backtrace(int sig, siginfo_t *info, void *dodgy)
{
  // Start by reporting on nature of signal
  std::string nature;
  if (sig == SIGSEGV) {
    nature = "Segmentation violation";
    switch (info->si_code) {
    case SEGV_MAPERR:
      nature += ", Address not mapped"; break;
    case SEGV_ACCERR:
      nature += ", Invalid permission for address"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Offending address = %p",info->si_addr);
    nature += tmp;
  }
  else if (sig == SIGBUS) {
    nature = "Bus error";
    switch (info->si_code) {
    case BUS_ADRALN:
      nature += ", Invalid address alignment"; break;
    case BUS_ADRERR:
      nature += ", Non-existent physical address"; break;
    case BUS_OBJERR:
      nature += ", Object specific hardware error"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Offending address = %p",info->si_addr);
    nature += tmp;
  }
  else if (sig == SIGFPE) {
    switch (info->si_code) {
    case FPE_INTDIV:
      nature += ", Integer divide-by-zero"; break;
    case FPE_INTOVF:
      nature += ", Integer overflow"; break;
    case FPE_FLTDIV:
      nature += ", Floating point divide-by-zero"; break;
    case FPE_FLTOVF:
      nature += ", Floating point overflow"; break;
    case FPE_FLTUND:
      nature += ", Floating point underflow"; break;
    case FPE_FLTRES:
      nature += ", Floating point loss of precision"; break;
    case FPE_FLTINV:
      nature += ", Invalid floating point operation"; break;
    case FPE_FLTSUB:
      nature += ", Subscript out of range"; break;
    default:
      nature += ", Unknown reason";
    }
    char tmp[128];
    sprintf(tmp,", Address of offending operation = %p",info->si_addr);
    nature += tmp;
  }
  std::cout << nature << std::endl;

  // Unwind stack
  void *stack_addr[MAX_SDEPTH];
  char **stack_strings;
  size_t stack_depth = backtrace(stack_addr,MAX_SDEPTH);
  stack_strings = backtrace_symbols(stack_addr,stack_depth);

  // Demangle C++ function names. Prior to that the strings from
  // backtrace_symbols have to be "cleaned up", a non-portable
  // and empirical abomination.  
  // The "cleaning up" of the strings below is based on empirical 
  // tests on jalapeno and jalapeno00 running RedHat 4.1.2
  for (unsigned int i=2; i<stack_depth; i++) {
    char bls[256]; strcpy(bls,stack_strings[i]);
    char *lp=NULL; char *lsp=NULL;
    for (lp=bls; *lp!='(' && *lp!='\0'; lp++) ; 
    if (*lp=='\0') {
      std::cout << std::string(bls) << std::endl;
    }
    else {
      *lp='\0'; lp++;
      for (lsp=lp; *lsp!='+'; lsp++) ; *lsp='\0';
      for (; *lsp!='['; lsp++) ;  
      size_t length; int status;
      char *demangled = abi::__cxa_demangle(lp,NULL,&length,&status);
      if (demangled) {
	std::cout << std::string(bls) << " " << std::string(demangled) << " " << std::string(lsp) << std::endl;
	free(demangled);
      }
      else { 
	std::cout << std::string(bls) << " " << std::string(lp) << " " << std::string(lsp) << std::endl;
      }
    }
  }
  
  free(stack_strings);

  exit(EXIT_FAILURE);
  return;
}
#endif // #ifdef __linux__

#endif // #ifndef stack_dump_h
