include(FetchContent)

# Silence the dev warning and make <PackageName>_ROOT variables work as intended.
cmake_policy(SET CMP0144 NEW)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Download prebuilt OR-Tools (C++ package) ----
#if(WIN32)
#  set(ORTOOLS_URL
#    "https://github.com/google/or-tools/releases/download/v9.12/or-tools_x64_VisualStudio2022_cpp_v9.12.4544.zip"
#  )
#elseif(UNIX AND NOT APPLE)
#  set(ORTOOLS_URL
#    "https://github.com/google/or-tools/releases/download/v9.12/or-tools_amd64_ubuntu-22.04_cpp_v9.12.4544.tar.gz"
#  )
#else()
#  message(FATAL_ERROR "Unsupported platform")
#endif()
set(ORTOOLS_VERSION "v9.15")
set(ORTOOLS_BUILD   "v9.15.6755")

if(WIN32)
  set(ORTOOLS_URL
    "https://github.com/google/or-tools/releases/download/${ORTOOLS_VERSION}/or-tools_x64_VisualStudio2022_cpp_${ORTOOLS_BUILD}.zip"
  )
elseif(UNIX AND NOT APPLE)
  set(ORTOOLS_URL
    "https://github.com/google/or-tools/releases/download/${ORTOOLS_VERSION}/or-tools_amd64_ubuntu-22.04_cpp_${ORTOOLS_BUILD}.tar.gz"
  )
else()
  message(FATAL_ERROR "Unsupported platform")
endif()

FetchContent_Declare(ortools_pkg URL ${ORTOOLS_URL})
FetchContent_MakeAvailable(ortools_pkg)

set(ORTOOLS_PREFIX "${ortools_pkg_SOURCE_DIR}")

# IMPORTANT:
# Make CMake search OR-Tools' prefix for *both* ortools and its dependencies (absl, protobuf, etc.)
list(PREPEND CMAKE_PREFIX_PATH "${ORTOOLS_PREFIX}" "${ORTOOLS_PREFIX}/lib/cmake")

# Optional: provide the modern <PackageName>_ROOT hint (this is what CMP0144 is about)
set(ortools_ROOT "${ORTOOLS_PREFIX}")

# ---- Load OR-Tools package config ----
find_package(ortools CONFIG REQUIRED)

# ---- Provide one stable target name for your project ----
add_library(cpsat::cpsat INTERFACE IMPORTED)
if(TARGET ortools::ortools)
  target_link_libraries(cpsat::cpsat INTERFACE ortools::ortools)
elseif(TARGET ortools)
  target_link_libraries(cpsat::cpsat INTERFACE ortools)
else()
  message(FATAL_ERROR "Found ortools package but no expected target (ortools::ortools or ortools).")
endif()
