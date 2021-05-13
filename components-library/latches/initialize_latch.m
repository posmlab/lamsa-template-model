%% function to initialize a new latch object

function latch = initialize_latch(app)

    name = app.LatchTabGroup.SelectedTab.Title;
    index = app.LatchParamIndices(name);
    vals = cell(1, index(2)-index(1)+1);
    for ii = index(1):index(2)
        val = str2double(app.LatchEditFields(ii).Value);
        vals{ii+1-index(1)} = val;
    end
    %latch = vals;
    latch = feval(name, vals{:});

end