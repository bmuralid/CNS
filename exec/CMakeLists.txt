add_executable("cns" "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_link_libraries("cns"
    PRIVATE
    ${EXEC_NAME})

