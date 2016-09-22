
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

set( general_sources 
  "helper_functions.f90" 
  "data_parameters.f90"  
  "eddington.f90"  
  "twostream.f90"
  "twostream_maxrand.f90"
  "nina_maxrand_wrapper.f90"
  "nina_maxrand.c"
  "gauss.c"
  )

add_library( maxrand STATIC ${general_sources} )
TARGET_LINK_LIBRARIES( maxrand ${EXTLIB} ${USERLIB} ${LIBS})

add_executable( "solve" "main.f90" )
target_link_libraries( "solve" maxrand )
