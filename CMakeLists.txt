cmake_minimum_required(VERSION 3.7.0)
include(CMakePrintHelpers)

# set the project name
project(FastJetContrib VERSION 1.101)
# set a pre-release version if relevant, e.g. "-beta1"
#set(PROJECT_VERSION_PRERELEASE "")
set(PROJECT_VERSION_PRERELEASE "-devel")

if(MSVC)
  execute_process(COMMAND 
                    git apply --whitespace=fix 0001-required-source-changes-for-windows-build.patch
                  WORKING_DIRECTORY
                    ${CMAKE_CURRENT_SOURCE_DIR}
                 )
endif()

#----------------------------------------------------------------------
# version information
set(FASTJETCONTRIB_VERSION "${PROJECT_VERSION}${PROJECT_VERSION_PRERELEASE}")

set(FASTJETCONTRIB_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(FASTJETCONTRIB_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(FASTJETCONTRIB_VERSION_PATCHLEVEL "${PROJECT_VERSION_PATCH}")
set(FASTJETCONTRIB_VERSION_PRERELEASE "${PROJECT_VERSION_PRERELEASE}")
# build a unique version number
if (FASTJETCONTRIB_VERSION_MINOR LESS 10)
  set(FASTJETCONTRIB_VERSION_MINOR_FOR_NUMBER "0${FASTJETCONTRIB_VERSION_MINOR}")
else()
  set(FASTJETCONTRIB_VERSION_MINOR_FOR_NUMBER "${FASTJETCONTRIB_VERSION_MINOR}")
endif()
if (FASTJETCONTRIB_VERSION_PATCHLEVEL LESS 10)
  set(FASTJETCONTRIB_VERSION_PATCHLEVEL_FOR_NUMBER "0${FASTJETCONTRIB_VERSION_PATCHLEVEL}")
else()
  set(FASTJETCONTRIB_VERSION_PATCHLEVEL_FOR_NUMBER "${FASTJETCONTRIB_VERSION_PATCHLEVEL}")
endif()
set(FASTJETCONTRIB_VERSION_NUMBER "${FASTJETCONTRIB_VERSION_MAJOR}${FASTJETCONTRIB_VERSION_MINOR_FOR_NUMBER}${FASTJETCONTRIB_VERSION_PATCHLEVEL_FOR_NUMBER}")

# print out the project name and version
cmake_print_variables(PROJECT_NAME FASTJETCONTRIB_VERSION)

# need this to be able to set the installation prefix in fastjet-config
include(GNUInstallDirs)

#----------------------------------------------------------------------
# basic C++ checks
#----------------------------------------------------------------------

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
if (MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()


if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM") # Intel compilers are used with -O3
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wshadow")
elseif (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:__cplusplus /utf-8 /W4 /Ox")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wshadow")
endif()

# all contribs require fastjet so get it now!
if(NOT TARGET fastjet::fastjet)
  find_package(fastjet REQUIRED)
endif()

# automatically determine the list of contribs, as done in configure
file(GLOB MAYBE_CONTRIB_DIRS LIST_DIRECTORIES true RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*")
foreach(MAYBE_CONTRIB_DIR ${MAYBE_CONTRIB_DIRS})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${MAYBE_CONTRIB_DIR} AND NOT ${MAYBE_CONTRIB_DIR} STREQUAL "Template")
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${MAYBE_CONTRIB_DIR}/VERSION
             OR EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${MAYBE_CONTRIB_DIR}/FJCONTRIB.cfg
          )
          list(APPEND fjcontrib_DIRNAMES ${MAYBE_CONTRIB_DIR})
        endif()
    endif()
endforeach()

foreach (FJCONTRIB_NAME ${fjcontrib_DIRNAMES})
  string(TOUPPER ${FJCONTRIB_NAME} FJCONTRIB_NAME_UPPER)
  option(FASTJETCONTRIB_ENABLE_${FJCONTRIB_NAME_UPPER} "Enable the ${FJCONTRIB_NAME} Contrib [default=ON]" ON)
endforeach()
option(FASTJETCONTRIB_ENABLE_FRAGILELIB "Enable fragile combined mono-library [default=OFF]" OFF)

set(fjcontrib_combined_TARGETS)

foreach (FJCONTRIB_NAME ${fjcontrib_DIRNAMES})

  string(TOUPPER ${FJCONTRIB_NAME} FJCONTRIB_NAME_UPPER)
  string(TOLOWER ${FJCONTRIB_NAME} FJCONTRIB_NAME_LOWER)
  if (FASTJETCONTRIB_ENABLE_${FJCONTRIB_NAME_UPPER})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/defines.h.in ${CMAKE_CURRENT_BINARY_DIR}/${FJCONTRIB_NAME}/include/${FJCONTRIB_NAME_LOWER}_defines.h)

    set(SRCLINE)
    set(INCLLINE)
    set(DEPLINE)
    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/Makefile FJCONTRIB_MAKEFILE)
    while(FJCONTRIB_MAKEFILE)
      list(POP_FRONT FJCONTRIB_MAKEFILE MAKEFILE_LINE)
      #message("${MAKEFILE_LINE}")
      if(MAKEFILE_LINE MATCHES "^SRCS=")
        string(REPLACE "SRCS=" "" SRCLINE "${MAKEFILE_LINE}")
        set(${FJCONTRIB_NAME_LOWER}_SOURCES "${SRCLINE}")
      endif()
      if(MAKEFILE_LINE MATCHES "^INSTALLED_HEADERS=")
        string(REPLACE "INSTALLED_HEADERS=" "" INCLLINE "${MAKEFILE_LINE}")
        set(${FJCONTRIB_NAME_LOWER}_INCLUDES "${INCLLINE}")
      endif()
      if(MAKEFILE_LINE MATCHES "^DEPENDS_ON\ =")
        string(REPLACE "DEPENDS_ON\ =" "" DEPLINE "${MAKEFILE_LINE}")
        string(STRIP "${DEPLINE}" DEPLINE)
        set(${FJCONTRIB_NAME_LOWER}_DEPENDENCIES "${DEPLINE}")
      endif()
    endwhile()

    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/incfastjet_CMakeLists.txt.in ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include/CMakeLists.txt COPYONLY)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/inccontrib_CMakeLists.txt.in ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include/fastjet/CMakeLists.txt COPYONLY)
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/incdir_CMakeLists.txt.in ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include/fastjet/contrib/CMakeLists.txt)
    endif()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/contrib_CMakeLists.txt.in ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/CMakeLists.txt)

    add_subdirectory(${FJCONTRIB_NAME})
    if (TARGET ${FJCONTRIB_NAME})
      if (FASTJETCONTRIB_ENABLE_FRAGILELIB)
        set_target_properties(${FJCONTRIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE True)
      endif()
      if (MSVC)
        target_compile_definitions(${FJCONTRIB_NAME} PUBLIC _USE_MATH_DEFINES)
        if (FASTJETCONTRIB_ENABLE_FRAGILELIB)
          target_compile_definitions(${FJCONTRIB_NAME} PRIVATE fastjetcontrib_EXPORTS)
        else()
          target_compile_definitions(${FJCONTRIB_NAME} PRIVATE ${FJCONTRIB_NAME_LOWER}_EXPORTS)
        endif()
      endif()
      list(APPEND fjcontrib_combined_TARGETS $<TARGET_OBJECTS:${FJCONTRIB_NAME}>)
    endif()    
  endif()

endforeach()

if (fjcontrib_combined_TARGETS AND FASTJETCONTRIB_ENABLE_FRAGILELIB)
   add_library(fastjetcontribfragile SHARED ${fjcontrib_combined_TARGETS})
   target_link_libraries(fastjetcontribfragile PUBLIC fastjet::fastjettools)
   foreach (FJCONTRIB_NAME ${fjcontrib_DIRNAMES})
     target_include_directories(fastjetcontribfragile
                                PUBLIC
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${FJCONTRIB_NAME}/include> 
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include>
                                  $<INSTALL_INTERFACE:include>
                                )
   endforeach()

   if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
     set_target_properties(fastjetcontribfragile PROPERTIES INSTALL_RPATH "@loader_path")
   elseif (NOT MSVC)
     set_target_properties(fastjetcontribfragile PROPERTIES INSTALL_RPATH "$ORIGIN")
   endif()

   add_library(fastjet::contrib::fastjetcontribfragile ALIAS fastjetcontribfragile)
   install(TARGETS
            fastjetcontribfragile
          EXPORT 
            FastjetContribTargets
          RUNTIME
            COMPONENT Runtime
            DESTINATION ${CMAKE_INSTALL_BINDIR}
          LIBRARY
            COMPONENT Runtime
            DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE
            COMPONENT Development
            DESTINATION ${CMAKE_INSTALL_LIBDIR}
         )
endif()

# allow Fastjet Contrib to work with find_package
export(EXPORT FastjetContribTargets
  NAMESPACE fastjet::contrib::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/fastjetcontribTargets.cmake"
)

install(EXPORT FastjetContribTargets
  NAMESPACE fastjet::contrib::
  FILE fastjetcontribTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastjetcontrib
)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/fastjetcontribConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/fastjetcontrib"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/fastjetcontribConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/fastjetcontribConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/fastjetcontribConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastjetcontrib
  )

