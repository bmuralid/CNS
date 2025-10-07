set (SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})

set (CNS_OBJECTS
    ${SRCDIR}/initEB.cpp
    ${SRCDIR}/cns.cpp
    ${SRCDIR}/AmrCoreCNS.cpp
    )

set (OBJECTS_SERIAL  ${CNS_OBJECTS} CACHE INTERNAL "Objects for serial build")
set (OBJECTS_MPI ${CNS_OBJECTS} CACHE INTERNAL "Objects for MPI build")

if(${USE_MPI})
   add_library(${EXEC_NAME} STATIC ${OBJECTS_MPI})
   set_target_properties (${EXEC_NAME} PROPERTIES COMPILE_DEFINITIONS "USE_MPI")
   target_link_libraries(${EXEC_NAME} PUBLIC AMReX::amrex)
   target_include_directories(${EXEC_NAME}
       PUBLIC
       ${MPI_INCLUDE_PATH}
       ${MPI_LIBRARY_PATH}
       ${CMAKE_CURRENT_SOURCE_DIR}/include
       ${CMAKE_CURRENT_SOURCE_DIR})
else()
   add_library(${EXEC_NAME} STATIC ${OBJECTS_MPI})
   target_link_libraries(${EXEC_NAME} PUBLIC AMReX::amrex)
   target_include_directories(${EXEC_NAME}
       PUBLIC
       ${CMAKE_CURRENT_SOURCE_DIR}/include
       ${CMAKE_CURRENT_SOURCE_DIR})
endif()

message(${EXEC_NAME})
