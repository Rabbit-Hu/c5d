include(CPM)

set(BOOST_INCLUDE_LIBRARIES "thread;asio;regex;filesystem;timer")

CPMAddPackage(
    NAME Boost
    URL "https://github.com/boostorg/boost/releases/download/boost-1.83.0/boost-1.83.0.tar.xz"
)
