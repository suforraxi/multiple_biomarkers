% apply a montage to all the strips present in the data (passed as a function)
% data      - fieldtrip data structure 
% fun       - handler to one of the possible montage function
%             the signature of the generic function takes 
%             two input parameters (channel label [cell], matrix of data [ch X samples])
%        
% outdata   - fieldtrip data structure with the data transformed according to the montage funtion fun 
function outdata = apply_montage2data_strip(data,fun)

outdata = data;
ntrial = numel(data.trial);


for i = 1: ntrial
    
   
    [stripslabel,dataXstrip]  = fun(data.label,data.trial{i});
    merged_data  = [];
    merged_label = [];
    for j = 1 : numel(dataXstrip)
         
        merged_data  = [merged_data ; dataXstrip{j} ];
        merged_label = [merged_label ; stripslabel{j}];   

    end
    outdata.trial{i} = merged_data;
    outdata.label    = merged_label;

    if(isempty(outdata.trial{i}))
        outdata = [];
        return
    end
end


    