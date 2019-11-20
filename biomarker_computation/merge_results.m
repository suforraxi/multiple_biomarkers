function outres = merge_results(outres_grid,outres_strip)

outres = [];
    
if (isempty(outres_grid) && ~isempty(outres_strip))
    
    outres = outres_strip;
    
     for i = 1 : numel(outres_grid.bio)       
            
           
            aux.strip       = outres_strip.extra{i};              
           
            outres.extra{i} = aux;
                    
    end
    
    
end

if(~isempty(outres_grid) && isempty(outres_strip))

    outres = outres_grid;
    
       for i = 1 : numel(outres_grid.bio)       
            
           
            aux.grid       = outres_grid.extra{i};              
           
            outres.extra{i} = aux;
                    
    end
    
end

if(~isempty(outres_grid) && ~isempty(outres_strip))
    
    outres        = outres_grid;
    outres.label  = {outres_grid.label{:} outres_strip.label{:}};   
    
    for i = 1 : numel(outres_grid.bio)       
            
            outres.bio{i}      =  [outres_grid.bio{i} ; outres_strip.bio{i}];
            
            aux.grid        = outres_grid.extra{i};
            aux.strip       = outres_strip.extra{i};              
            outres.extra{i} = aux;
                    
    end

end