function WriteLinopToFile(FN,OpsCell)
if(~iscell(OpsCell))
    OpsCell={OpsCell};
end
fid=fopen(FN,'wt');
for i=1:numel(OpsCell)
    if(iscell(OpsCell{i}))
        for j=1:numel(OpsCell{i})
            fprintf(fid,"%s\n",OpsCell{i}{j});
        end
        if(i<numel(OpsCell))
            fprintf(fid,"nextlinop\n");
        end
    else
        fprintf(fid,"%s\n",OpsCell{i});
    end
end
fclose(fid);