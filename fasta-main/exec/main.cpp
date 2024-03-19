#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <fasta/Fasta.h>
#include <fasta/Matrix.h>
#include <fasta/version.h>
#include <fasta/utils/Timers.h>
#include <fasta/utils/Counters.h>

// #if defined(__linux__) || defined(__APPLE)
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
// #endif


auto main(int argc, char** argv) -> int {

using namespace std;
using namespace timer;
using namespace counter;
using namespace fasta;

  // #if defined(__linux__) || defined(__APPLE)
  signal(SIGSEGV, handler);   // install our handler
  // #endif

  try {

    cxxopts::Options options(*argv,
      "FASTA - an exact Nearest Neighbor Search algorithm for dense data using cosine similarity.");

    std::string qmat_fpath;
    std::string dmat_fpath;
    std::string output_fpath;
    std::string method;
    std::string mode;
    std::string scale_mode;

    // set up searcher config options
    fasta::Config config;

    // clang-format off
    options
      .add_options()
      ("h,help", "Show help")
      ("v,version", "Print the current version number")
      ("q,query", "Query matrix file to find neighbors for.", cxxopts::value(qmat_fpath), "FILE")
      ("d,database", "Database matrix file to find neighbors in.", cxxopts::value(dmat_fpath), "FILE")
      ("o,output", "Output file", cxxopts::value(output_fpath), "FILE")
      ("m,method", "Method to use for searching (gemm, l2, bsi, or fasta).",
        cxxopts::value(method)->default_value("gemm"), "S")
      ("s,mode", "Search mode (sim or knn).", cxxopts::value(mode)->default_value("sim"), "S")
      ("e,epsilon", "Minimum search similarity", cxxopts::value(config.epsilon)->default_value("0.9"), "F")
      ("k,nnbrs", "Maximum number of neighbors", cxxopts::value(config.max_k)->default_value("10"), "N")
      ("t,nthreads", "Number of threads", cxxopts::value(config.nthreads)->default_value("4"), "N")
      ("n,normalize", "Normalize input matrices", cxxopts::value<bool>()->default_value("false"))
      ("p,permute", "Permute columns in L2 methods", cxxopts::value<bool>()->default_value("false"))
      ("scale", "Scale input matrices", cxxopts::value(scale_mode)->default_value(""))
      ("io", "I/O mode: read query matrix and write it to output file",
        cxxopts::value<bool>()->default_value("false"))
      ("info", "Info mode: read query and/or database matrices and provide some basic info about them.",
        cxxopts::value<bool>()->default_value("false"))
      ("debug", "Enable debugging", cxxopts::value<bool>()->default_value("false"))
      ("l,logpath", "Path to log file.", cxxopts::value<std::string>()->default_value(""))
    ;
    // clang-format on

    auto opt = options.parse(argc, argv);

    if (opt.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
    }

    if (opt.count("version")) {
      std::cout << "fasta, version " << version.major << "." << version.minor << "." << version.patch << std::endl;
      return 0;
    }

    // set up logging
    auto logger = spdlog::stderr_color_mt("stdout");   // register standard logger to stdout
    auto logpath = opt["logpath"].as<std::string>();
    if (!logpath.empty()) {
      try
      {
        // register rotating file logger to given log file path
        auto max_size = 1048576 * 5;
        auto max_files = 3;
        auto file_logger = spdlog::rotating_logger_mt("file_logger", logpath, max_size, max_files);
      }
      catch (const spdlog::spdlog_ex &ex)
      {
        std::cerr << "Log init failed: " << ex.what() << std::endl;
      }
      spdlog::info("Rotating file logger added at {}.", logpath);
    }
    spdlog::flush_every(std::chrono::seconds(3)); // periodically flush all *registered* loggers every 3 seconds:
    spdlog::set_level(spdlog::level::info); // Set global log level to info
    if (opt.count("debug")) {
      spdlog::set_level(spdlog::level::debug); // Set global log level to debug
      spdlog::enable_backtrace(32); // Store the latest 32 messages in a buffer. Older messages will be dropped.
    }
    Timer timers[NTIMERS];
    Counter counters[NCOUNTERS];

    // info mode
    if(opt.count("info")){
      timers[EXEC].start();
      if(!qmat_fpath.empty()){
        spdlog::info("Read query matrix {}.", qmat_fpath);
        auto qmat = Matrix<float>::read(qmat_fpath, 'd');
        cout << qmat->get_info() << endl;
        delete qmat;
      }
      if(!dmat_fpath.empty()){
        spdlog::info("Read database matrix {}.", qmat_fpath);
        auto dmat = Matrix<float>::read(dmat_fpath, 'd');
        cout << dmat->get_info() << endl;
        delete dmat;
      }
      timers[EXEC].stop();
      cout << "Execution time: " << timers[EXEC] << endl;
      return 0;
    }

    // I/O mode
    if(qmat_fpath.empty()){
      std::cerr << "The query file is required." << std::endl << options.help() << std::endl;
      return 1;
    }
    if(opt.count("io")){
      if(output_fpath.empty()){
        std::cerr << "The output file is required in I/O mode." << std::endl << options.help() << std::endl;
        return 1;
      }
      timers[EXEC].start();
      spdlog::info("Read query matrix {}.", qmat_fpath);
      auto qmat = Matrix<float>::read(qmat_fpath, 'd');
      if(!scale_mode.empty()){
        spdlog::info("Scaling query matrix using {} scaling mode.", scale_mode);
        qmat->scale(scale_mode);
      }
      if(opt.count("normalize")){
        spdlog::info("Normalizing query matrix.");
        qmat->normalize('2');
      }
      spdlog::info("Writing query matrix to {}.", output_fpath);
      qmat->write(output_fpath);
      delete qmat;
      timers[EXEC].stop();
      cout << "Execution time: " << timers[EXEC] << endl;
      return 0;
    }

    if(dmat_fpath.empty()){
      std::cerr << "The database file is required." << std::endl << options.help() << std::endl;
      return 1;
    }

    // ensure valid method and mode -- TODO: move to utils
    auto mode_id = NMODES;
    auto method_id = NMETHODS;
    transform(method.begin(), method.end(), method.begin(), ::tolower);
    transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
    for(uint8_t i=0; i < NMETHODS; ++i){
      if(MethodNames[i] == method){
        method_id = Methods(i);
      }
    }
    if(method_id == NMETHODS){
      std::cerr << "Invalid method." << std::endl << options.help() << std::endl;
      return 1;
    }
    for(uint8_t i=0; i < NMODES; ++i){
      if(ModeNames[i] == mode){
        mode_id = Modes(i);
      }
    }
    if(mode_id == NMODES){
      std::cerr << "Invalid search mode." << std::endl << options.help() << std::endl;
      return 1;
    }
    config.mode = mode_id;
    config.method = method_id;
    config.permute_columns = opt.count("permute");
    spdlog::info("Search method: {}", MethodNames[method_id]);
    spdlog::info("Search mode: {}", ModeNames[mode_id]);
    spdlog::info("Permute columns: {}", opt.count("permute") ? "on" : "off");

    timers[EXEC].start();
    spdlog::debug("Execution timer started");
    // read input matrices
    timers[READ].start();
    spdlog::debug("Read timer started");
    auto qmat = Matrix<float>::read(qmat_fpath, 'd');
    auto dmat = Matrix<float>::read(dmat_fpath, 'd');
    timers[READ].stop();
    spdlog::debug("Read complete");
    spdlog::info("Database: " + dmat->get_info());
    spdlog::info("Query: " + qmat->get_info());

    // pre-process data
    if(!scale_mode.empty()){
      spdlog::info("Scaling database matrix using {} scaling mode.", scale_mode);
      double * factors;
      dmat->scale(scale_mode, &factors);
      // apply db scaling factors to qmat
      spdlog::info("Applying database matrix scaling factors to query matrix.");
      qmat->scale(scale_mode, factors);
      free(factors);
    }

    // set up searcher
    timers[PREPROCESS].start();
    Fasta<float> searcher(config, *dmat);
    timers[PREPROCESS].stop();
    // execute search
    if(mode_id == KNN){
      auto nnbrs = opt["nnbrs"].as<idx_t>();
      Matrix<sim_t> res(qmat->nrows, nnbrs);
      timers[SEARCH].start();
      searcher.knn_search(*qmat, nnbrs, res);
      timers[SEARCH].stop();
      // write results
      if(!output_fpath.empty()){
        timers[WRITE].start();
        res.write(output_fpath);
        timers[WRITE].stop();
      }
    } else {
      csr_t res;
      timers[SEARCH].start();
      searcher.search(*qmat, res);
      timers[SEARCH].stop();
      // write results
      if(!output_fpath.empty()){
        timers[WRITE].start();
        res.write(output_fpath);
        timers[WRITE].stop();
      }
    }

    // transfer index timers and counters
    timers[EXEC].stop();
    add_timers(searcher.get_index()->timers, timers);
    add_counters(searcher.get_index()->counters, counters);
    cout << endl << "Counts:" << endl << get_counts_string(counters, "\t");
    cout << endl << "Overall Times:" << endl << get_times_string(timers, "\t");
    delete qmat;
    delete dmat;
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
  return 0;
}
