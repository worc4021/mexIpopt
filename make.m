LDFLAGS = '-bundle -Wl,-twolevel_namespace -undefined error -stdlib=libc++';
BINROOT = '/Users/Manuel/Documents/Development/binaries';


mex('-output','ipopt',...
    ['-I',fullfile(BINROOT,'include','coin-or')],...
    ['-L',fullfile(BINROOT,'lib')],...
    ['CXXFLAGS=-std=c++17'],...
    ['LDFLAGS=',LDFLAGS],...
    ['CC=/usr/bin/xcrun  gcc'],...
    ['CXX=/usr/bin/xcrun  g++'],...
    ['LD=/usr/bin/xcrun  gcc'],...
    ['LDXX=/usr/bin/xcrun  g++'],...
    '-lipopt',...
    'main.cpp')