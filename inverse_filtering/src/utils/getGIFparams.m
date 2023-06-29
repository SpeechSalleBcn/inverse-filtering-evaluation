function GIFparams = getGIFparams(C)

% Add analysis information to the table
T = table;
T.FileName = C.filesData.speechFileName;
fileInfo = C.filesData.fileInfo;
fileFields = C.filesData.fileFields;

for i=1:length(fileFields)
    if ~isnan(str2double(fileInfo.(fileFields{i})))
        T.(fileFields{i}) = str2double(fileInfo.(fileFields{i}));
    elseif ~isempty(fileInfo.(fileFields{i}))
        T.(fileFields{i}) = fileInfo.(fileFields{i});
    else
        T.(fileFields{i}) = "";
    end
end

% TODO: afegir columna de rang d'f0 (f0r) a nivell de codi matlab o ho fem a l'script de python??
% if (T.f0<C.gif.f0bins(2))
%     f0r = C.gif.f0label{1};
% elseif (T.f0<C.gif.f0bins(3))
%     f0r = C.gif.f0label{2};
% else
%     f0r = C.gif.f0label{3};
% end
%T.f0
%f0r
% verificar que T.f0 es correspon amb f0r

fileRowTable = C.paramsTable(strcmp(C.paramsTable.GIFmethod, string(C.gif.GIFMethod)),:);
for i=1:length(fileFields)
    if ~ismissing(string(eval(['fileRowTable.', fileFields{i}])))
        fileRowTable = fileRowTable(strcmp( string(eval(['fileRowTable.', fileFields{i}])) , string(eval(['T.', fileFields{i}]))), :);
    end
end

GIFparams.VTorder = fileRowTable.VTorder;
GIFparams.GSorder = fileRowTable.GSorder;
GIFparams.LipRad = fileRowTable.LipRad;
if strcmp(C.gif.method,'Original-IAIF') || strcmp(C.gif.method,'IOP-IAIF')
    GIFparams.HPflag = fileRowTable.HPflag;
elseif strcmp(C.gif.method,'QCP')
    GIFparams.DQ = fileRowTable.DQ;
    GIFparams.PQ = fileRowTable.PQ;
    GIFparams.RQ = fileRowTable.RQ;
    GIFparams.STcompensation = fileRowTable.STcompensation;
end
end