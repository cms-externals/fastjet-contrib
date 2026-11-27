cmake_minimum_required(VERSION 3.7.0)
include(CMakePrintHelpers)

# set the project name
project(FastJetContrib VERSION 1.105)
# set a pre-release version if relevant, e.g. "-beta1"
set(PROJECT_VERSION_PRERELEASE "")
#set(PROJECT_VERSION_PRERELEASE "-devel")

# conda-forge applies this patch as a recipe source patch and does not ship the file
if(MSVC AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/0001-required-source-changes-for-windows-build.patch)
  # keep git from treating an enclosing repository as this tree's, where 'git apply' skips every file and exits 0
  get_filename_component(FJCONTRIB_WINPATCH_CEILING ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
  set(FJCONTRIB_GIT ${CMAKE_COMMAND} -E env GIT_CEILING_DIRECTORIES=${FJCONTRIB_WINPATCH_CEILING} git)
  execute_process(COMMAND
                    ${FJCONTRIB_GIT} apply --whitespace=fix 0001-required-source-changes-for-windows-build.patch
                  WORKING_DIRECTORY
                    ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE FJCONTRIB_WINPATCH_RESULT
                  ERROR_VARIABLE FJCONTRIB_WINPATCH_ERROR
                 )
  if(NOT FJCONTRIB_WINPATCH_RESULT EQUAL 0)
    # a re-configure finds the patch already applied
    execute_process(COMMAND
                      ${FJCONTRIB_GIT} apply --reverse --check 0001-required-source-changes-for-windows-build.patch
                    WORKING_DIRECTORY
                      ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE FJCONTRIB_WINPATCH_APPLIED
                    OUTPUT_QUIET ERROR_QUIET
                   )
    if(NOT FJCONTRIB_WINPATCH_APPLIED EQUAL 0)
      message(FATAL_ERROR "0001-required-source-changes-for-windows-build.patch does not apply:\n${FJCONTRIB_WINPATCH_ERROR}")
    endif()
    message(STATUS "fjcontrib: windows source patch already applied")
  else()
    message(STATUS "fjcontrib: applied windows source patch")
  endif()
endif()

#----------------------------------------------------------------------
# version information
set(FASTJETCONTRIB_VERSION "${PROJECT_VERSION}${PROJECT_VERSION_PRERELEASE}")

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
file(GLOB fjcontrib_DIRNAMES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     "${CMAKE_CURRENT_SOURCE_DIR}/*/VERSION" "${CMAKE_CURRENT_SOURCE_DIR}/*/FJCONTRIB.cfg")
string(REGEX REPLACE "/[^;]*" "" fjcontrib_DIRNAMES "${fjcontrib_DIRNAMES}")
list(REMOVE_DUPLICATES fjcontrib_DIRNAMES)
list(REMOVE_ITEM fjcontrib_DIRNAMES Template)

option(FASTJETCONTRIB_ENABLE_FRAGILELIB "Enable fragile combined mono-library [default=OFF]" OFF)

set(fjcontrib_combined_TARGETS)
set(fjcontrib_INSTALL_TARGETS)

foreach (FJCONTRIB_NAME ${fjcontrib_DIRNAMES})

  string(TOUPPER ${FJCONTRIB_NAME} FJCONTRIB_NAME_UPPER)
  string(TOLOWER ${FJCONTRIB_NAME} FJCONTRIB_NAME_LOWER)
  option(FASTJETCONTRIB_ENABLE_${FJCONTRIB_NAME_UPPER} "Enable the ${FJCONTRIB_NAME} Contrib [default=ON]" ON)
  if (NOT FASTJETCONTRIB_ENABLE_${FJCONTRIB_NAME_UPPER})
    continue()
  endif()
  set(FJCONTRIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME})
  set(FJCONTRIB_BIN ${CMAKE_CURRENT_BINARY_DIR}/${FJCONTRIB_NAME})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/cmake-templates/defines.h.in ${FJCONTRIB_BIN}/include/${FJCONTRIB_NAME_LOWER}_defines.h)

  # sources, headers and dependencies come from the contrib's own Makefile
  set(FJCONTRIB_SRCS)
  set(FJCONTRIB_INSTALLED_HEADERS)
  set(FJCONTRIB_DEPENDS_ON)
  file(STRINGS ${FJCONTRIB_SRC}/Makefile FJCONTRIB_MAKEFILE REGEX "^(SRCS|INSTALLED_HEADERS|DEPENDS_ON)[ \t]*=")
  foreach(MAKEFILE_LINE ${FJCONTRIB_MAKEFILE})
    string(REGEX MATCH "^([A-Z_]+)[ \t]*=(.*)" MAKEFILE_LINE "${MAKEFILE_LINE}")
    string(REGEX MATCHALL "[^ \t]+" FJCONTRIB_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
  endforeach()
  string(REGEX REPLACE "[^;]+" "${FJCONTRIB_SRC}/\\0" FJCONTRIB_SRCS "${FJCONTRIB_SRCS}")

  if (EXISTS ${FJCONTRIB_SRC}/include)
    string(REGEX REPLACE "[^;]+" "${FJCONTRIB_SRC}/include/fastjet/contrib/\\0" FJCONTRIB_INSTALLED_HEADERS "${FJCONTRIB_INSTALLED_HEADERS}")
  else()
    string(REGEX REPLACE "[^;]+" "${FJCONTRIB_SRC}/\\0" FJCONTRIB_INSTALLED_HEADERS "${FJCONTRIB_INSTALLED_HEADERS}")
    file(COPY ${FJCONTRIB_INSTALLED_HEADERS} DESTINATION ${FJCONTRIB_BIN}/include/fastjet/contrib)
  endif()
  install(FILES
            ${FJCONTRIB_INSTALLED_HEADERS}
            ${FJCONTRIB_BIN}/include/${FJCONTRIB_NAME_LOWER}_defines.h
          DESTINATION
            include/fastjet/contrib
         )

  if (NOT FJCONTRIB_SRCS)
    continue()
  endif()
  if (FASTJETCONTRIB_ENABLE_FRAGILELIB)
    add_library(${FJCONTRIB_NAME} OBJECT ${FJCONTRIB_SRCS})
    set_target_properties(${FJCONTRIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE True)
  else()
    add_library(${FJCONTRIB_NAME} SHARED ${FJCONTRIB_SRCS})
    # per-contrib output directories keep the build-tree layout of one directory per contrib
    set_target_properties(${FJCONTRIB_NAME} PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${FJCONTRIB_BIN}
                          RUNTIME_OUTPUT_DIRECTORY ${FJCONTRIB_BIN}
                          ARCHIVE_OUTPUT_DIRECTORY ${FJCONTRIB_BIN})
  endif()
  target_link_libraries(${FJCONTRIB_NAME} PUBLIC fastjet::fastjettools ${FJCONTRIB_DEPENDS_ON})
  if (FJCONTRIB_NAME STREQUAL "SDFPlugin")
    target_link_libraries(${FJCONTRIB_NAME} PUBLIC fastjet::fastjetplugins)
  endif()
  target_include_directories(${FJCONTRIB_NAME}
                             PUBLIC
                              $<BUILD_INTERFACE:${FJCONTRIB_SRC}/include>
                              $<BUILD_INTERFACE:${FJCONTRIB_BIN}/include>
                              $<INSTALL_INTERFACE:include>
                             )

  if (MSVC)
    target_compile_definitions(${FJCONTRIB_NAME} PUBLIC _USE_MATH_DEFINES)
    if (FASTJETCONTRIB_ENABLE_FRAGILELIB)
      target_compile_definitions(${FJCONTRIB_NAME} PRIVATE fastjetcontrib_EXPORTS)
    else()
      target_compile_definitions(${FJCONTRIB_NAME} PRIVATE ${FJCONTRIB_NAME_LOWER}_EXPORTS)
    endif()
  endif()

  if (FASTJETCONTRIB_ENABLE_FRAGILELIB)
    list(APPEND fjcontrib_combined_TARGETS $<TARGET_OBJECTS:${FJCONTRIB_NAME}>)
  else()
    add_library(fastjet::contrib::${FJCONTRIB_NAME} ALIAS ${FJCONTRIB_NAME})
    list(APPEND fjcontrib_INSTALL_TARGETS ${FJCONTRIB_NAME})
  endif()

endforeach()

if (fjcontrib_combined_TARGETS AND FASTJETCONTRIB_ENABLE_FRAGILELIB)
   add_library(fastjetcontribfragile SHARED ${fjcontrib_combined_TARGETS})
   target_link_libraries(fastjetcontribfragile PUBLIC fastjet::fastjettools)
   # SDFPlugin's headers construct fastjet plugins inline, as its per-contrib target links
   if (TARGET SDFPlugin)
     target_link_libraries(fastjetcontribfragile PUBLIC fastjet::fastjetplugins)
   endif()
   foreach (FJCONTRIB_NAME ${fjcontrib_DIRNAMES})
     target_include_directories(fastjetcontribfragile
                                PUBLIC
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${FJCONTRIB_NAME}/include> 
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${FJCONTRIB_NAME}/include>
                                )
   endforeach()
   target_include_directories(fastjetcontribfragile PUBLIC $<INSTALL_INTERFACE:include>)

   add_library(fastjet::contrib::fastjetcontribfragile ALIAS fastjetcontribfragile)
   list(APPEND fjcontrib_INSTALL_TARGETS fastjetcontribfragile)
endif()

# header-only contribs install no library, and an export set needs at least one
if (fjcontrib_INSTALL_TARGETS)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${fjcontrib_INSTALL_TARGETS} PROPERTIES INSTALL_RPATH "@loader_path")
  elseif (NOT MSVC)
    set_target_properties(${fjcontrib_INSTALL_TARGETS} PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()
  install(TARGETS
            ${fjcontrib_INSTALL_TARGETS}
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
endif()

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

