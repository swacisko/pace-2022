# ============================================================
# GoogleTest + Google Benchmark dependencies (WSL2-friendly)
# ============================================================

include(FetchContent)

# ---------------- GoogleTest ----------------
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(gtest_force_shared_crt OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)

# ---------------- Google Benchmark ----------------
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.8.5
)

# Download & make available
FetchContent_MakeAvailable(googletest benchmark)

# Enable CTest + GTest helpers
enable_testing()
include(GoogleTest)

# Exposed targets:
#   GTest::gtest
#   GTest::gtest_main
#   benchmark::benchmark
#   benchmark::benchmark_main