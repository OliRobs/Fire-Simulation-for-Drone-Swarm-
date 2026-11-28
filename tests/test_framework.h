#pragma once
// ============================================================================
// Minimal test framework for the headless tests (no third-party dependency).
//
// TEST(name) { ... } defines a test and registers it at start-up; CHECK(cond) records a failure and carries
// on, so one run reports every failing check. SLOW_TEST(name) is the same for a test that takes most of a
// second or more: --fast leaves those out (ctest gives them their own entry, labelled `slow`).
// The tests live in one file per module (tests/test_*.cpp) and are all linked into one executable, wf_tests;
// tests/test_main.cpp runs them and takes filters:
//
//   wf_tests                    run every test
//   wf_tests <substring>        run the tests whose name contains <substring>
//   wf_tests --file <substring> run the tests whose source file name contains <substring> (ctest uses this)
//   wf_tests --fast             leave out the SLOW_TESTs (combine with the filters above)
//   wf_tests --slow             run only the SLOW_TESTs
//   wf_tests --list             print every test's name and file, run nothing
// ============================================================================
#include <cstdio>
#include <vector>

namespace wf_test {
    // One registered test: its name, the source file it is in (__FILE__), its function and whether it is slow
    struct TestCase {
        const char* name;
        const char* file;
        void (*fn)();
        bool slow;
    };

    // Every test, in registration order (registration order between files is up to the linker;
    // test_main.cpp sorts by file so the run order is stable)
    std::vector<TestCase>& registry();

    // Failed CHECKs so far in this run
    int& failures();

    // Registers a test from a static object's constructor (see TEST)
    struct Registrar {
        Registrar(const char* name, const char* file, void (*fn)(), bool slow) {
            registry().push_back({ name, file, fn, slow });
        }
    };
}

// Defines `static void name()` and registers it under the name "name"
#define TEST(name)                                                                  \
    static void name();                                                             \
    static ::wf_test::Registrar reg_##name(#name, __FILE__, name, false);           \
    static void name()

// As TEST, for a test that takes most of a second or more (left out by --fast)
#define SLOW_TEST(name)                                                             \
    static void name();                                                             \
    static ::wf_test::Registrar reg_##name(#name, __FILE__, name, true);            \
    static void name()

// Records a failure (with file and line) when cond is false, and carries on with the test
#define CHECK(cond) do { if (!(cond)) { ++::wf_test::failures();                    \
    std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)
