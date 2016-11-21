set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")

# sets HAVE_FMEMOPEN
add_subdirectory(ext)

add_library (ffindex ffindex.c ffutil.c)

target_include_directories (ffindex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (ffindex_shared SHARED ffindex.c ffutil.c)

if(NOT HAVE_FMEMOPEN)
	target_link_libraries(ffindex ext)
	target_link_libraries(ffindex_shared ext)
endif()

add_executable(ffindex_apply
	ffindex_apply.c
)
target_link_libraries (ffindex_apply ffindex)

add_executable(ffindex_apply_mpi
    ffindex_apply_mpi.c
)

target_link_libraries (ffindex_apply_mpi ffindex)

if(HAVE_MPI)
find_package(MPI)
if(MPI_C_FOUND)
    add_definitions("-DHAVE_MPI=1")
    add_subdirectory(mpq)
	set_property(TARGET ffindex_apply_mpi PROPERTY COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
	set_property(TARGET ffindex_apply_mpi PROPERTY LINK_FLAGS "${MPI_C_LINK_FLAGS}")
    target_include_directories(ffindex_apply_mpi PUBLIC ${MPI_C_INCLUDE_PATH})
    target_link_libraries (ffindex_apply_mpi mpq ${MPI_C_LIBRARIES})
endif()
endif()

install(PROGRAMS
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_apply_mpi
    DESTINATION bin
)



add_executable(ffindex_build
  ffindex_build.c
)
target_link_libraries (ffindex_build ffindex)


add_executable(ffindex_from_fasta
  ffindex_from_fasta.c
)
target_link_libraries (ffindex_from_fasta ffindex)


add_executable(ffindex_get
  ffindex_get.c
)
target_link_libraries (ffindex_get ffindex)


add_executable(ffindex_modify
  ffindex_modify.c
)
target_link_libraries (ffindex_modify ffindex)


add_executable(ffindex_unpack
  ffindex_unpack.c
)
target_link_libraries (ffindex_unpack ffindex)


add_executable(ffindex_order
  ffindex_order.c
)
target_link_libraries (ffindex_order ffindex)


add_executable(ffindex_from_fasta_with_split
    ffindex_from_fasta_with_split.c
)
target_link_libraries (ffindex_from_fasta_with_split ffindex)


install(PROGRAMS 
	ffindex.h 
	ffutil.h
	DESTINATION include
)

install(TARGETS ffindex
  LIBRARY DESTINATION lib64 COMPONENT libraries
	ARCHIVE DESTINATION lib64 COMPONENT libraries
)


install(TARGETS ffindex_shared
  LIBRARY DESTINATION lib64 COMPONENT libraries
	ARCHIVE DESTINATION lib64 COMPONENT libraries
)


INSTALL(PROGRAMS 
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_apply
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_build
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_from_fasta
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_get
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_modify
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_unpack
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_order
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ffindex_from_fasta_with_split
  DESTINATION bin
)
