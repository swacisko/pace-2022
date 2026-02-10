include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- OR-Tools build options (minimal) ----
set(BUILD_DEPS          ON  CACHE BOOL "" FORCE)
set(BUILD_SAMPLES       OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES      OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING       OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON        OFF CACHE BOOL "" FORCE)
set(BUILD_JAVA          OFF CACHE BOOL "" FORCE)
set(BUILD_DOTNET        OFF CACHE BOOL "" FORCE)
set(BUILD_CXX           ON  CACHE BOOL "" FORCE)

# ---- Download OR-Tools ----
FetchContent_Declare(
  ortools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        v9.10
)

FetchContent_MakeAvailable(ortools)

## ---- Example CP-SAT target (optional) ----
#add_executable(cpsat_demo main.cpp)
#target_link_libraries(cpsat_demo PRIVATE ortools)
