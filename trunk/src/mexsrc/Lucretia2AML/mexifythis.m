% Compile Lucretia2AML mex file

% Set Paths UAP and XERCES:
% (you also need to set your LD_LIBRARY_PATH to include $XERCESROOT/lib)
% XERCESROOT='/nfs/slac/g/nlcrd/u01/smolloy/xerces-c-src_2_8_0';
% if strcmpi(getenv('PROC_TYPE'),'x86_64')
%   UAPROOT='/a/sulky18/u02/smolloy/accelerator-ml/uap/trunk64';
% else
%   UAPROOT='/a/sulky18/u02/smolloy/accelerator-ml/uap/trunk';
% end
XERCESROOT='/home/smolloy/Documents/xerces-c-src_2_8_0';
UAPROOT='/home/smolloy/accelerator-ml/uap/trunk';

mex('-v', ...
  ['-L',UAPROOT,'/lib -lantlr -luap'], ...
  ['-L',XERCESROOT,'/lib -lxerces-c'], ...
  ['-I',XERCESROOT,'/include'], ...
  ['-I',UAPROOT,'/ANTLR270/antlr'], ...
  ['-I',UAPROOT,'/ANTLR270'], ...
  ['-I',UAPROOT], ...
  ['-I',UAPROOT,'/AML'], ...
  ['-I',UAPROOT,'/UAP'], ...
  ['-I',UAPROOT,'/Bmad'], ...
  'Lucretia2AML.cpp', 'LEle2UAPNode.cpp', 'make_amlmarker.cpp', 'make_amlsben.cpp', 'make_amlquad.cpp', 'addName.cpp', ...
  'addS.cpp', 'addL.cpp', 'addB.cpp', 'addBerr.cpp', 'LEle2UAPNode_ol.cpp', 'unsplitquad.cpp', 'make_amlaper.cpp', ...
  'make_amlorient.cpp', 'unsplitsben.cpp', 'make_amldrift.cpp', 'make_amlbpm.cpp', 'make_amlxcor.cpp', ...
  'addmover.cpp', 'unsplitxcor.cpp', 'unsplitycor.cpp', 'unsplitbpm.cpp', ...
  'CreateAMLController.cpp', 'make_amlycor.cpp', 'make_amlsext.cpp', 'CreateAMLGirder.cpp', 'unsplitsext.cpp', ...
  'addFloodLand.cpp', 'equalNodes.cpp', 'make_amloct.cpp',...
  [UAPROOT,'/ANTLR270/Parser.o'], ...
  [UAPROOT,'/ANTLR270/CharScanner.o'], ...
  [UAPROOT,'/ANTLR270/TreeParser.o'], ...
  [UAPROOT,'/ANTLR270/LLkParser.o'], ...
  [UAPROOT,'/ANTLR270/NoViableAltForCharException.o'], ...
  [UAPROOT,'/ANTLR270/NoViableAltException.o']);

clear XERCESROOT UAPROOT

