% Run this script after running test.m to animate the resulting solution.
% If not using test.m, use solve_lamsa_se (or similar) with syntax 
% [sol, t_times] = solve_lamsa_se(tspan*, loading_motor, unlatching_motor*, load, latch*, spring).
% Variable names with asterisks can be changed as you like, but this script
% assumes the others will be named as above.
% Running the script will play an animation and save it to the variable
% "frames"; it can be replayed by calling movie(frames) or, for better
% control, by using the VideoViewer app. 

v = VideoWriter('test');
open(v);
% Set-up stuff
L4 = spring.rest_length + loading_motor.rest_length;
L3 = load.lengths(3);
L2 = load.lengths(2);
LVec = [L2, L3, L4];
l0 = spring.rest_length + loading_motor.rest_length;
theta0 = load.theta_0;

times = sol(:, 1);
UL_time = t_times(1); % unlatching time

% Interpolate times and angles so timesteps between frames are constant
thetas = sol(:,2);
ints = linspace(times(1), times(end), 150);
angle = interp1(times, thetas, ints, 'spline');

rawY1s = sol(:,12);
y1s = interp1(times, rawY1s, ints, 'spline');

rawSs = sol(:, 4);
Ss = interp1(times, rawSs, ints, 'spline');

frames = struct('cdata', cell(1, length(ints)),'colormap',cell(1, length(ints)));

for k = 1:length(angle)
    % Lever arm
    qLoad = quiver(-L2*cos(angle(k)), -L2*sin(angle(k)), (L2+L3)*cos(angle(k)), (L2+L3)*sin(angle(k)), '-');
    qLoad.ShowArrowHead = 0;
    hold on;

    % Vectors for spring/muscle stuff
    fVec = [L2*(cos(theta0) - cos(angle(k))), l0 + L2*(sin(theta0) - sin(angle(k)))];
    uVec = fVec/norm(fVec);
    mVec = (loading_motor.rest_length - y1s(k))*uVec;
    sVec = fVec - mVec;

    %%%% ------- Uncomment this and comment out qMuscle and qSpring
    %%%% portions to ignore muscle/spring distinction; be warned that this
    %%%% will hide any weird behaviour with y1 ------------- %%%%%%%%%%%%%%%%%%
%     qSpringMuscle = quiver(-L2*cos(theta0), -L2*sin(theta0) - l0, fVec(1), fVec(2), 'LineStyle', '-', 'Color', 'k', 'AutoScale','off');
%     qSpringMuscle.ShowArrowHead = 0;
    %%%% --------------------------------------------%%%%%%%%%%%%%%%%%%%

    qMuscle = quiver(-L2*cos(theta0), -L2*sin(theta0) - l0, mVec(1), mVec(2), 'LineStyle', '-', 'Color', 'r', 'AutoScale','off');
    qMuscle.ShowArrowHead = 0;
    
    qSpring = quiver(-L2*cos(theta0) + mVec(1), -L2*sin(theta0) - l0 + mVec(2), sVec(1), sVec(2), 'LineStyle', '-', 'Color', 'g', 'AutoScale','off');
    qSpring.ShowArrowHead = 0;

    plot(-L2*cos(theta0), -L2*sin(theta0) - l0, 'k.', 'MarkerSize', 6);
    plot(0, 0, 'k.', 'MarkerSize', 6);


%     plot(-L2*cos(theta0) - Ss(k), -L2*sin(theta0), 'k.', 'MarkerSize', 6);  % Adds crude latch animation
   
    % Writes out time and latched status
    txt = "Time: " + num2str(ints(k), 4);
    txt = txt + newline;
    if ints(k) > UL_time
        txt = txt + "Unlatched";
    else
        txt = txt + "Latched";
    end
    text(0.8*max(L2, L3), max(L2, L3), txt);

    axis.style = 'square';

    set(gca, 'XLim', [-1.5*max(LVec), 1.5*max(LVec)], 'YLim', [-1.5*max(LVec), 1.5*max(LVec)]);

    frames(k) = getframe;
    writeVideo(v,frames(k));
    hold off;
    drawnow
end
close(v);

% Convert to gif
% !ffmpeg -i test.avi -vf "fps=30,palettegen=stats_mode=diff" palette.png
% !ffmpeg -i test.avi -i palette.png -lavfi "fps=10,paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" test.gif