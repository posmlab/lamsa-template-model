%% function to initialize a new component

function obj = initialize_component(name, param_indices, editfields)

    index = param_indices(name);
    vals = cell(1, index(2)-index(1)+1);
    for ii = index(1):index(2)
        val = str2double(editfields(ii).Value);
        vals{ii+1-index(1)} = val;
    end
    obj = feval(name, vals{:});

end