BINROOT = '/Users/manuel/Documents/Development/binaries';
MEXBUILDINGBLOCKS = '/Users/manuel/Documents/Development/MexBuildingBlocks';

LDIR = '"/MATLAB Drive/build/lib"';
IDIR = '"/MATLAB Drive/build/include/coin-or"';
MEXBUILDINGBLOCKS = fullfile(pwd,'../MexBuildingBlocks');



if ismac
    mex('-output','ipopt',...
        ['-I',fullfile(BINROOT,'include','coin-or')],...
        ['-L',fullfile(BINROOT,'lib')],...
        ['-I',MEXBUILDINGBLOCKS],...
        '-lipopt',...
        'main.cpp')
elseif isunix

    mex('-output','ipopt','-v',...
        '-R2018a',...
        ['-I',IDIR],...
        ['-L',LDIR],...
        ['-I',MEXBUILDINGBLOCKS],...
        '"/MATLAB Drive/build/lib/libipopt.a"',...
        '"/MATLAB Drive/build/lib/libcoinhsl.a"',...
        '-lmwblas',...
        '-lmwlapack',...
        '-lgfortran',...
        'main.cpp')
end
