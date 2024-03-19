/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once


#include <exception>
#include <string>
#include <utility>
#include <vector>
#include <sstream>

#ifdef __GNUG__
#include <cxxabi.h>
#endif

namespace fasta {


/// Base class for Fasta exceptions
class FastaException : public std::exception {
   public:
    explicit FastaException(const std::string& msg);

    FastaException(
            const std::string& msg,
            const char* funcName,
            const char* file,
            int line);

    /// from std::exception
    const char* what() const noexcept override;

    std::string msg;
};

FastaException::FastaException(const std::string& m) : msg(m) {}

FastaException::FastaException(
        const std::string& m,
        const char* funcName,
        const char* file,
        int line) {
    int size = snprintf(
            nullptr,
            0,
            "Error in %s at %s:%d: %s",
            funcName,
            file,
            line,
            m.c_str());
    msg.resize(size + 1);
    snprintf(
            &msg[0],
            msg.size(),
            "Error in %s at %s:%d: %s",
            funcName,
            file,
            line,
            m.c_str());
}

const char* FastaException::what() const noexcept {
    return msg.c_str();
}


/// Handle multiple exceptions from worker threads, throwing an appropriate
/// exception that aggregates the information
/// The pair int is the thread that generated the exception
void handleExceptions(
        std::vector<std::pair<int, std::exception_ptr>>& exceptions) {
    if (exceptions.size() == 1) {
        // throw the single received exception directly
        std::rethrow_exception(exceptions.front().second);

    } else if (exceptions.size() > 1) {
        // multiple exceptions; aggregate them and return a single exception
        std::stringstream ss;

        for (auto& p : exceptions) {
            try {
                std::rethrow_exception(p.second);
            } catch (std::exception& ex) {
                if (ex.what()) {
                    // exception message available
                    ss << "Exception thrown from index " << p.first << ": "
                       << ex.what() << "\n";
                } else {
                    // No message available
                    ss << "Unknown exception thrown from index " << p.first
                       << "\n";
                }
            } catch (...) {
                ss << "Unknown exception thrown from index " << p.first << "\n";
            }
        }

        throw FastaException(ss.str());
    }
}

/** bare-bones unique_ptr
 * this one deletes with delete [] */
template <class T>
struct UniqueScope {
    const T* ptr;
    explicit UniqueScope(const T* ptr = nullptr) : ptr(ptr) {}
    void release() {
        ptr = nullptr;
    }
    void set(const T* ptr_in) {
        ptr = ptr_in;
    }
    void swap(UniqueScope<T>& other) {
        std::swap(ptr, other.ptr);
    }
    ~UniqueScope() {
        delete[] ptr;
    }
};

/** same but deletes with the simple delete (least common case) */
template <class T>
struct UniqueScopeSingle {
    const T* ptr;
    explicit UniqueScopeSingle(const T* ptr = nullptr) : ptr(ptr) {}
    void release() {
        ptr = nullptr;
    }
    void set(const T* ptr_in) {
        ptr = ptr_in;
    }
    void swap(UniqueScopeSingle<T>& other) {
        std::swap(ptr, other.ptr);
    }
    ~UniqueScopeSingle() {
        delete ptr;
    }
};

// From
// https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
// Make typeids more readable
std::string demangle_cpp_symbol(const char* name) {
#ifdef __GNUG__
    int status = -1;
    const char* res = abi::__cxa_demangle(name, nullptr, nullptr, &status);
    std::string sres;
    if (status == 0) {
        sres = res;
    }
    free((void*)res);
    return sres;
#else
    // don't know how to do this on other platforms
    return std::string(name);
#endif
}

} // namespace fasta
