#include <csignal>
#include <cstdarg>
#include <iostream>
#include <sstream>

#include <Omega_h_fail.hpp>

extern "C" void Omega_h_signal_handler(int s);

namespace Omega_h {

#ifdef OMEGA_H_THROW
exception::exception(std::string const& msg_in) : msg(msg_in) {}

const char* exception::what() const noexcept { return msg.c_str(); }
#endif

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

void fail(char const* format, ...) {
  va_list vlist;
  va_start(vlist, format);
#ifdef OMEGA_H_THROW
  char buffer[2048];
  std::vsnprintf(buffer, sizeof(buffer), format, vlist);
  va_end(vlist);
  throw Omega_h::exception(buffer);
#else
  std::vfprintf(stderr, format, vlist);
  va_end(vlist);
  std::abort();
#endif
}

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic pop
#endif

static struct {
  int code;
  const char* name;
} const known_signals[] = {
#ifndef _MSC_VER
    {SIGSYS, "bad system call"},
    {SIGTSTP, "terminal stop"},
    {SIGQUIT, "quit"},
    {SIGHUP, "hangup"},
#endif
    {SIGABRT, "abort"}, {SIGTERM, "termination"},
    {SIGSEGV, "segmentation fault"}, {SIGINT, "interrupt"},
    {SIGILL, "illegal instruction"}, {SIGFPE, "floating point exception"}};
constexpr std::size_t NSIGS =
    (sizeof(known_signals) / sizeof(known_signals[0]));

void protect() {
  for (std::size_t i = 0; i < NSIGS; ++i) {
    signal(known_signals[i].code, Omega_h_signal_handler);
  }
}
}  // namespace Omega_h

extern "C" void Omega_h_signal_handler(int s) {
  static volatile sig_atomic_t already_dying = 0;
  if (already_dying) return;
  already_dying = 1;
  std::stringstream ss;
  for (std::size_t i = 0; i < Omega_h::NSIGS; ++i) {
    if (s == Omega_h::known_signals[i].code) {
      ss << "Omega_h caught signal: " << Omega_h::known_signals[i].name << "\n";
    }
  }
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
  signal(s, SIG_DFL);
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang diagnostic pop
#endif
  ::raise(s);
}
