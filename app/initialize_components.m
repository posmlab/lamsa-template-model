
%% function to initialize all components

function [spring, loading_motor, unlatching_motor, load, latch] = initialize_components(app)

    spring_name = app.SpringTabGroup.SelectedTab.Title;
    latch_name = app.LatchTabGroup.SelectedTab.Title;
    mass_name = app.MassTabGroup.SelectedTab.Title;
    loading_motor_name = app.LoadingMotorTabGroup.SelectedTab.Title;
    unlatching_motor_name = app.UnlatchingMotorTabGroup.SelectedTab.Title;

    spring = initialize_component(spring_name, app.SpringParamIndices, app.SpringEditFields);
    loading_motor = initialize_component(loading_motor_name, app.LoadingMotorParamIndices, app.LoadingMotorEditFields);
    unlatching_motor = initialize_component(unlatching_motor_name, app.UnlatchingMotorParamIndices, app.UnlatchingMotorEditFields);
    load = initialize_component(mass_name, app.MassParamIndices, app.MassEditFields);
    latch = initialize_component(latch_name, app.LatchParamIndices, app.LatchEditFields);

end

%% helper function to initialize a single new component

function obj = initialize_component(name, param_indices, editfields)

    index = param_indices(name);
    vals = cell(1, index(2)-index(1)+1);
    for ii = index(1):index(2)
        val = editfields(ii).Value;
        vals{ii+1-index(1)} = val;
    end
    obj = feval(name, vals{:});

end