% Make list of all shared object libs for compiled mex binaries and copy
% into Lucretia Libs directory
exclude={'libstdc++'};
if ismac
  [~,solibs1]=system(sprintf('find ../* -name "*.%s" | xargs otool -L | grep "/Users"',mexext));
  [~,solibs2]=system(sprintf('find ../* -name "*.%s" | xargs otool -L | grep "/usr/local"',mexext));
  solibs=[solibs1 solibs2];
else
  [~,solibs]=system(sprintf('find ../* -name "*.%s" | xargs ldd | grep "=> /"',mexext));
end
solibs=unique(regexp(strsplit(solibs,'\n'),'/(\S+)','match','once'));
solibs=solibs(cellfun(@(x) ~isempty(x),solibs));
for ilib=1:length(solibs)
  fprintf('Copying %s to ../Libs/%s ...\n',solibs{ilib},mexext)
  try
    copyfile(solibs{ilib},sprintf('../Libs/%s',mexext),'f');
    system(sprintf('svn add ../Libs/%s/%s',mexext,regexprep(solibs{ilib},'.*/','')));
  catch
  end
end
