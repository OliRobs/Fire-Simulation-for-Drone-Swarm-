// ============================================================================
// Test runner for wf_tests: runs the registered tests (test_framework.h), optionally filtered.
//
//   wf_tests                    every test
//   wf_tests <substring>        tests whose name contains <substring>
//   wf_tests --file <substring> tests whose source file name contains <substring>
//   wf_tests --fast             leave out the SLOW_TESTs;  --slow: run only them
//   wf_tests --list             list the tests (name and file) without running them
//
// Exit code: 0 if every check passed, 1 if any failed, 2 for a bad command line or a name/file filter that
// matches no test (so a mistyped ctest entry fails instead of passing with nothing run). A file whose tests
// are all slow, run with --fast, passes with nothing run.
// ============================================================================
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include "test_framework.h"

namespace wf_test {
    // Function-local statics: safe to use from the Registrar constructors in any file, whatever the order
    // in which the files' static objects are built
    std::vector<TestCase>& registry() {
        static std::vector<TestCase> r;
        return r;
    }
    int& failures() {
        static int f = 0;
        return f;
    }
}

namespace {
    // The file name without its directories, for --list and --file (__FILE__ may be a full path)
    std::string baseName(const char* path) {
        const char* slash = std::max(std::strrchr(path, '/'), std::strrchr(path, '\\'));
        return slash ? slash + 1 : path;
    }
}

int main(int argc, char** argv) {
    using namespace wf_test;

    // --- Command line ---
    bool list = false, fastOnly = false, slowOnly = false;
    std::string nameFilter, fileFilter;
    for (int a = 1; a < argc; ++a) {
        if (std::strcmp(argv[a], "--list") == 0) list = true;
        else if (std::strcmp(argv[a], "--fast") == 0) fastOnly = true;
        else if (std::strcmp(argv[a], "--slow") == 0) slowOnly = true;
        else if (std::strcmp(argv[a], "--file") == 0 && a + 1 < argc) fileFilter = argv[++a];
        else if (argv[a][0] != '-' && nameFilter.empty()) nameFilter = argv[a];
        else {
            std::fprintf(stderr, "usage: wf_tests [<name substring>] [--file <file substring>] [--fast | --slow] "
                                 "[--list]\n");
            return 2;
        }
    }

    // --- Selection: stable order, files alphabetically, then registration order within each file ---
    std::vector<TestCase> tests = registry();
    std::stable_sort(tests.begin(), tests.end(),
                     [](const TestCase& a, const TestCase& b) { return baseName(a.file) < baseName(b.file); });
    tests.erase(std::remove_if(tests.begin(), tests.end(), [&](const TestCase& t) {
                    return (!nameFilter.empty() && std::strstr(t.name, nameFilter.c_str()) == nullptr)
                        || (!fileFilter.empty() && baseName(t.file).find(fileFilter) == std::string::npos);
                }), tests.end());
    if (tests.empty()) {
        std::fprintf(stderr, "No test matches the filter\n");
        return 2;
    }
    // Speed filter after the name/file check, so "--file x --fast" on a file of slow tests is not an error
    tests.erase(std::remove_if(tests.begin(), tests.end(), [&](const TestCase& t) {
                    return (fastOnly && t.slow) || (slowOnly && !t.slow);
                }), tests.end());
    if (tests.empty()) {
        std::printf("No test left after --fast / --slow\n");
        return 0;
    }

    if (list) {
        for (const TestCase& t : tests)
            std::printf("%s  (%s%s)\n", t.name, baseName(t.file).c_str(), t.slow ? ", slow" : "");
        std::printf("%zu test(s)\n", tests.size());
        return 0;
    }

    // --- Run: name, any output and FAIL lines from the test, then ok/FAILED with the time taken ---
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    for (const TestCase& t : tests) {
        const int before = failures();
        std::printf("%s\n", t.name);
        const auto t0 = Clock::now();
        t.fn();
        const double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        std::printf("  %s (%.0f ms)\n", failures() == before ? "ok" : "FAILED", ms);
    }
    const double total = std::chrono::duration<double>(Clock::now() - start).count();
    std::printf("\n%s: %zu test(s), %d failure(s), %.1f s\n", failures() ? "FAILED" : "PASSED", tests.size(),
                failures(), total);
    return failures() ? 1 : 0;
}
