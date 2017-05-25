INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_LTE_CHANNELS lte_channels)

FIND_PATH(
    LTE_CHANNELS_INCLUDE_DIRS
    NAMES lte_channels/api.h
    HINTS $ENV{LTE_CHANNELS_DIR}/include
        ${PC_LTE_CHANNELS_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    LTE_CHANNELS_LIBRARIES
    NAMES gnuradio-lte_channels
    HINTS $ENV{LTE_CHANNELS_DIR}/lib
        ${PC_LTE_CHANNELS_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LTE_CHANNELS DEFAULT_MSG LTE_CHANNELS_LIBRARIES LTE_CHANNELS_INCLUDE_DIRS)
MARK_AS_ADVANCED(LTE_CHANNELS_LIBRARIES LTE_CHANNELS_INCLUDE_DIRS)

