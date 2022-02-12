
LDIR = '/home/manuel/build/lib';
IDIR = '/home/manuel/build/include/coin-or';
MEXBUILDINGBLOCKS = '/mnt/c/repos/MexBuildingBlocks';



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
        '/home/manuel/build/lib/libipopt.a',...
        '/home/manuel/build/lib/libcoinhsl.a',...
        '-lmwblas',...
        '-lmwlapack',...
        '-lgfortran',...
        'main.cpp')
end
