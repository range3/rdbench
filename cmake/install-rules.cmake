install(
    TARGETS rdbench_exe
    RUNTIME COMPONENT rdbench_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
