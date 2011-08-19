%
% script for mex-ing XSIF tools XSIFParse and GetXSIFDictionary
%

%==============================================================

%
% first get the OS environment variable
%
ostype = getenv('OSTYPE') ;
%
% do certain things differently in Solaris and Linux
%
if (strcmp(ostype,'solaris'))      % solaris
     xsifpath = '/afs/slac/g/nlc/codes/matliar/bin/xsif' ;
     xsiflib = 'xsif_sun' ;
elseif (strcmp(ostype,'linux'))
     xsifpath = '/afs/slac/g/nlc/codes/xsif/binlinux' ;
     xsiflib = 'xsif_linux' ;
end
%
% build the command strings
%
     mexXSIFParse = ['mex -I. -L.',...
			     ' -l',xsiflib,' XSIFParse.f'] ;
     mexXSIFDict = ['mex -I. -L.',...
			     ' -l',xsiflib,' GetXSIFDictionary.f'] ;
%
% execute the command strings
%
     eval(mexXSIFParse) ;
     eval(mexXSIFDict) ;
%
