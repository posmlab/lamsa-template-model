function showKinematics(axesHandle,evnt, xAxisVariable, yAxisVariable, app)

% ptH = getappdata(axesHandle,'CurrentPoint')
% pt = get(gcf,'CurrentPoint') %Getting click position

% properties(axesHandle)
xClick = evnt.IntersectionPoint(1);
yClick = evnt.IntersectionPoint(2);

%% using the user's click to calculate and use the closest pixel value. 
if strcmp(app.x_log_space.Value,'log')
    xrange = [log10(app.xmin.Value) log10(app.xmax.Value)];
    looping_values_x = logspace(xrange(1),xrange(2),app.n.Value);
    x_values_on_graph = arrayfun(@(x) log10(x), looping_values_x);
    clickDistances = arrayfun(@(listValue) abs(listValue - xClick), x_values_on_graph);
else
    xrange = [app.xmin.Value app.xmax.Value];
    looping_values_x = linspace(xrange(1),xrange(2),app.n.Value);
    clickDistances = arrayfun(@(listValue) abs(listValue - xClick), looping_values_x);
end
[minimum_distance_x, indexOfValue_x] = min(clickDistances);
xUsed = looping_values_x(indexOfValue_x);

[~, xindex] = ismember(xAxisVariable, app.IV1DropDown.Items);
xefIndex = ddIndexToefIndex(app, xindex, xAxisVariable);
changeEditField(app, xAxisVariable, xefIndex, xUsed);

if strcmp(app.y_log_space.Value,'log')
    yrange = [log10(app.ymin.Value) log10(app.ymax.Value)];
    looping_values_y = logspace(yrange(1),yrange(2),app.n.Value);
    y_values_on_graph = arrayfun(@(y) log10(y), looping_values_y);
    clickDistances = arrayfun(@(listValue) abs(listValue - yClick), y_values_on_graph);
    
else
    yrange = [app.ymin.Value app.ymax.Value];
    looping_values_y = linspace(yrange(1),yrange(2),app.n.Value);
    clickDistances = arrayfun(@(listValue) abs(listValue - yClick), looping_values_y);
end

[minimum_distance_y, indexOfValue_y] = min(clickDistances);
yUsed = looping_values_y(indexOfValue_y);

[~, yindex] = ismember(yAxisVariable, app.IV1DropDown.Items);
yefIndex = ddIndexToefIndex(app, yindex, yAxisVariable);
changeEditField(app, yAxisVariable, yefIndex, yUsed);

%% initializing LaMSA component structs
[spring, loading_motor, unlatching_motor, load, latch] = initialize_components(app);

% calling solve model
[sol,transition_times] = solve_model(loading_motor,unlatching_motor,load,latch,spring);

plotNames = {'Time vs. y coordinate of load', 'Time vs. y velocity of load', 'Time vs. y-forces on load',...
             'Time vs. x coordinate of latch','Time vs. x velocity of latch','Time vs. x-forces on latch'};

[indexs, booleanValue] = listdlg('ListString', plotNames, ...
                             'Name','Plot Selection',...
                             'PromptString','Which kinematic plot(s) would you like to see?',... 
                             'ListSize', [800, 500]);

if (booleanValue)
figure;

if ((ismember(1,indexs) || ismember(2,indexs) || ismember(3,indexs)) && ...
    (ismember(4,indexs) || ismember(5,indexs) || ismember(6,indexs)))
    numColumns = 2;  
else
    numColumns = 1;
end

numYPlots = ismember(1, indexs) + ismember(2, indexs) + ismember(3, indexs);
numXPlots = ismember(4, indexs) + ismember(5, indexs) + ismember(6, indexs);
numRows = max(numXPlots, numYPlots);


% I'm sorry. This is awful. 
% This is textbook how-NOT-to do if statements.
% but this is what seems readable and how my brain structured it in my
% head.

if (numColumns == 1)
    if (ismember(1,indexs) || ismember(2, indexs) || ismember(3, indexs)) % if plotting y column
        if (ismember(1, indexs))
            subplot(numRows, numColumns, 1)
            plot(sol(:,1),sol(:,2))
            xlabel("time",'Interpreter','latex')
            ylabel("y coordinate of load",'Interpreter','latex')
            title("y coordinate of load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(2, indexs))
            if (ismember(1,indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            plot(sol(:,1),sol(:,3))
            xlabel("time",'Interpreter','latex')
            ylabel("y velocity of load",'Interpreter','latex')
            title("y velocity of load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(3, indexs))
            if (ismember(1, indexs) && ismember(2, indexs))
                subplot(numRows, numColumns, 3)
            elseif (ismember(1, indexs) || ismember(2, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            y_forces = sol(:,10)+sol(:,9)-sol(:,7);
            plot(sol(:,1),y_forces)
            xlabel("time",'Interpreter','latex')
            ylabel("y forces on load",'Interpreter','latex')
            title("y forces on load vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        
    else % else: plotting x column
        if (ismember(4, indexs))
            subplot(numRows, numColumns, 1)
            plot(sol(:,1),sol(:,4))
            xlabel("time",'Interpreter','latex')
            ylabel("x coordinate of latch",'Interpreter','latex')
            title("x coordinate of latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(5,indexs))
            if (ismember(4, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            plot(sol(:,1),sol(:,5))
            xlabel("time",'Interpreter','latex')
            ylabel("x velocity of latch",'Interpreter','latex')
            title("x velocity of latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
        if (ismember(6, indexs))
            if (ismember(4, indexs) && ismember(5, indexs))
                subplot(numRows, numColumns, 3)
            elseif (ismember(4, indexs) || ismember(5, indexs))
                subplot(numRows, numColumns, 2)
            else
                subplot(numRows, numColumns, 1)
            end
            x_forces = sol(:,11)+sol(:,8)+sol(:,6);
            plot(sol(:,1),x_forces)
            xlabel("time",'Interpreter','latex')
            ylabel("x forces on latch",'Interpreter','latex')
            title("x forces on latch vs. time",'Interpreter','latex')
            set(gca,'TickLabelInterpreter','Latex');
        end
    end
end

if (numColumns == 2)
    if (ismember(1, indexs)) % t vs y
        subplot(numRows, numColumns,1)
        plot(sol(:,1),sol(:,2))
        xlabel("time",'Interpreter','latex')
        ylabel("y coordinate of load",'Interpreter','latex')
        title("y coordinate of load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(4, indexs)) % t vs x
        subplot(numRows, numColumns,2)
        plot(sol(:,1),sol(:,4))
        xlabel("time",'Interpreter','latex')
        ylabel("x coordinate of latch",'Interpreter','latex')
        title("x coordinate of latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(2, indexs)) % t vs ydot
        if (ismember(1, indexs))
            subplot(numRows, numColumns, 3)
        else
            subplot(numRows, numColumns, 1)
        end
        plot(sol(:,1),sol(:,3))
        xlabel("time",'Interpreter','latex')
        ylabel("y velocity of load",'Interpreter','latex')
        title("y velocity of load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
    if (ismember(5, indexs)) % t vs xdot
        if (ismember(4, indexs))
            subplot(numRows, numColumns, 4)
        else
            subplot(numRows, numColumns, 2)
        end
        plot(sol(:,1),sol(:,5))
        xlabel("time",'Interpreter','latex')
        ylabel("x velocity of latch",'Interpreter','latex')
        title("x velocity of latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end

    if (ismember(3, indexs)) % t vs yforces
        if ((~ismember(1, indexs)) && (~ismember(2, indexs)))
            subplot(numRows, numColumns, 1)
        elseif xor(ismember(1, indexs),ismember(2, indexs))
            subplot(numRows, numColumns, 3)
        else
            subplot(numRows, numColumns, 5)
        end
        y_forces = sol(:,10)+sol(:,9)-sol(:,7);
        plot(sol(:,1), y_forces)
        xlabel("time",'Interpreter','latex')
        ylabel("y forces on load",'Interpreter','latex')
        title("y forces on load vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end

    if (ismember(6, indexs)) % t vs xforces
        if ((~ismember(4, indexs)) && (~ismember(5, indexs)))
            subplot(numRows, numColumns, 2)
        elseif xor(ismember(4, indexs),ismember(5, indexs))
            subplot(numRows, numColumns, 4)
        else
            subplot(numRows, numColumns, 6)
        end
        x_forces = sol(:,11)+sol(:,8)+sol(:,6);
        plot(sol(:,1), x_forces)
        xlabel("time",'Interpreter','latex')
        ylabel("x forces on latch",'Interpreter','latex')
        title("x forces on latch vs. time",'Interpreter','latex')
        set(gca,'TickLabelInterpreter','Latex');
    end
end

end

end
