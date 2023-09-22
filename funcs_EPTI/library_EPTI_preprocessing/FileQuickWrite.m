function FileQuickWrite(FileName,X)
% write matrix to file, .txt 

fid = fopen(FileName,'wt');
ncol = size(X,2);
fprintf(fid, [repmat('%f ',1,ncol),'\r\n'], X');
fclose(fid);

end