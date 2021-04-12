cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
set(CMAKE_INSTALL_PREFIX "~" CACHE PATH "Install prefix" FORCE)




    SET(GSL_CONFIG_PREFER_PATH "$ENV{GSL_HOME}/bin" CACHE STRING "preferred path to GSL (gsl-config)")
    FIND_PROGRAM(GSL_CONFIG gsl-config
      ${GSL_CONFIG_PREFER_PATH}
      $ENV{LIB_DIR}/bin
      /usr/local/bin/
      /usr/bin/
      )
    # MESSAGE("DBG GSL_CONFIG ${GSL_CONFIG}")
    
    IF (GSL_CONFIG) 
      # set CXXFLAGS to be fed into CXX_FLAGS by the user:
      SET(GSL_CXX_FLAGS "`${GSL_CONFIG} --cflags`")
      
      # set INCLUDE_DIRS to prefix+include
      EXEC_PROGRAM(${GSL_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE GSL_PREFIX)
      SET(GSL_INCLUDE_DIR ${GSL_PREFIX}/include CACHE STRING INTERNAL)

      # set link libraries and link flags
      EXEC_PROGRAM(${GSL_CONFIG}
          ARGS --libs
          OUTPUT_VARIABLE GSL_LIBRARIES)
      
      ## extract link dirs for rpath  
      EXEC_PROGRAM(${GSL_CONFIG}
        ARGS --libs
        OUTPUT_VARIABLE GSL_CONFIG_LIBS )

      ## split off the link dirs (for rpath)
      ## use regular expression to match wildcard equivalent "-L*<endchar>"
      ## with <endchar> is a space or a semicolon
      STRING(REGEX MATCHALL "[-][L]([^ ;])+" 
        GSL_LINK_DIRECTORIES_WITH_PREFIX 
        "${GSL_CONFIG_LIBS}" )
        #      MESSAGE("DBG  GSL_LINK_DIRECTORIES_WITH_PREFIX=${GSL_LINK_DIRECTORIES_WITH_PREFIX}")

      ## remove prefix -L because we need the pure directory for LINK_DIRECTORIES
      
      IF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
        STRING(REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES ${GSL_LINK_DIRECTORIES_WITH_PREFIX} )
      ENDIF (GSL_LINK_DIRECTORIES_WITH_PREFIX)
      SET(GSL_EXE_LINKER_FLAGS "-Wl,-rpath,${GSL_LINK_DIRECTORIES}" CACHE STRING INTERNAL)
      #      MESSAGE("DBG  GSL_LINK_DIRECTORIES=${GSL_LINK_DIRECTORIES}")
      #      MESSAGE("DBG  GSL_EXE_LINKER_FLAGS=${GSL_EXE_LINKER_FLAGS}")

      #      ADD_DEFINITIONS("-DHAVE_GSL")
      #      SET(GSL_DEFINITIONS "-DHAVE_GSL")
      MARK_AS_ADVANCED(
        GSL_CXX_FLAGS
        GSL_INCLUDE_DIR
        GSL_LIBRARIES
        GSL_LINK_DIRECTORIES
        GSL_DEFINITIONS
      )
      
    ELSE(GSL_CONFIG)

      IF (GSL_FIND_REQUIRED)
         MESSAGE(FATAL_ERROR "Could not find gsl-config. Please set it manually. GSL_CONFIG=${GSL_CONFIG}")
      ELSE (GSL_FIND_REQUIRED)
         MESSAGE(STATUS "Could not find GSL")
         # TODO: Avoid cmake complaints if GSL is not found
      ENDIF (GSL_FIND_REQUIRED)

    ENDIF(GSL_CONFIG)





IF(GSL_LIBRARIES)
  IF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)

    SET(GSL_FOUND 1)
    
    MESSAGE(STATUS "Using GSL from ${GSL_PREFIX}")

  ENDIF(GSL_INCLUDE_DIR OR GSL_CXX_FLAGS)
ENDIF(GSL_LIBRARIES)



find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS} ./)
add_definitions(${ROOT_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")

set(ROOT_LIBRARIES ${ROOT_LIBRARIES})            
add_executable(daqana_v6f4 daqana_v6f4.cc RootFileManager.cc RootFileManager.hh)
target_link_libraries(daqana_v6f4 ${ROOT_LIBRARIES} fftw3 ${GSL_LIBRARIES} -lsnappy )
install (TARGETS daqana_v6f4 DESTINATION bin)

