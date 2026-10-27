#ifndef RCPPTRUST_R_CALLBACK_ERROR_H
#define RCPPTRUST_R_CALLBACK_ERROR_H

// Captures an error raised inside an R callback so it can be re-raised
// once the thread-safe solver core has unwound. The error is stored by
// value and re-thrown as a fresh exception rather than through
// std::exception_ptr/std::rethrow_exception: with clang + libc++ (CRAN's
// fedora-clang) the rethrown object no longer matched Rcpp's catch
// clauses and surfaced as "c++ exception (unknown reason)".
#include <Rcpp.h>
#include <exception>
#include <string>

class RCallbackError {
 public:
  // Call only from inside a catch block.
  void capture() {
    try {
      throw;
    } catch (Rcpp::internal::InterruptedException &) {
      kind_ = INTERRUPT;
    } catch (Rcpp::LongjumpException &e) {
      // The unwind token resumes the original R condition unchanged.
      kind_ = LONGJUMP;
      token_ = e.token;
    } catch (std::exception &e) {
      kind_ = MESSAGE;
      msg_ = e.what();
    } catch (...) {
      kind_ = MESSAGE;
      msg_ = "unknown C++ exception in R callback";
    }
  }

  bool failed() const { return kind_ != NONE; }

  // Re-raise the captured error (no-op if nothing was captured). Under
  // Rcpp's unwind protection an R-level interrupt arrives as a
  // LongjumpException, and callers only rethrow after failed(), so the
  // INTERRUPT and NONE cases are defensive and not reachable from tests.
  void rethrow() const {
    switch (kind_) {
      case INTERRUPT: throw Rcpp::internal::InterruptedException();  // # nocov
      case LONGJUMP: throw Rcpp::LongjumpException(token_);
      case MESSAGE: throw Rcpp::exception(msg_.c_str(), false);
      case NONE: break;  // # nocov start
    }
  }  // # nocov end

 private:
  enum Kind { NONE, INTERRUPT, LONGJUMP, MESSAGE };
  Kind kind_ = NONE;
  Rcpp::RObject token_;  // keeps the unwind token protected
  std::string msg_;
};

#endif
