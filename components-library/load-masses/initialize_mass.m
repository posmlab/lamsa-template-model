%% function to initialize a new mass object

function mass = initialize_mass(app)

    name = app.MassTabGroup.SelectedTab.Title;
    index = app.MassParamIndices(name);
    vals = cell(1, index(2)-index(1)+1);
    for ii = index(1):index(2)
        val = str2double(app.MassEditFields(ii).Value);
        vals{ii+1-index(1)} = val;
    end
    mass = feval(name, vals{:});

end