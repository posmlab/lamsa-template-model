function writeInfoToFile(m_eff, transition_times, sol, loading_motor,unlatching_motor,load,latch,spring, outputDirectory)

params = struct();

params.m_L = latch.mass;

params.m_eff = m_eff;

params.v0 = latch.v_0;

func_struct = functions(loading_motor.Force);
params.loading_motor = rmfield(func_struct,{'file','type','within_file_path'});

func_struct = functions(spring.Force);
params.spring = rmfield(func_struct,{'file','type','within_file_path'});

func_struct = functions(unlatching_motor.Force);
params.unlatching_force = rmfield(func_struct,{'file','type','within_file_path'});

func_struct = functions(latch.y_L{1});
params.latch_shape = rmfield(func_struct,{'file','type','within_file_path'});

params.unlatch_time = transition_times(1);

params.launch_time = transition_times(2);

%% Making output directory
if ~isfolder(outputDirectory)%checks for and possibly creates output directory
    mkdir(outputDirectory)
end

timeStampString = string(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy-MM-dd_HH-mm-ss-SSS'));

%% Writing .json output
pretty = prettyjson(['parameters: ' jsonencode(params)]);%makes the parameter json file readable

fileName = outputDirectory + sprintf("/parameters_%s.json", timeStampString);%ensures the file is in the output directory

fileID = fopen(fileName, 'w');

fwrite(fileID, pretty, 'char');
fclose(fileID);

%% save solution data to csv and json files
csvFilePath = outputDirectory + "/raw_data--" + timeStampString + ".csv";%same as above for json file, ensures location

headers = ["Time", "y", "ydot", "x", "xdot", "normal force on latch x", ...
    "normal force on load y", "frictional force on latch x", ...
    "frictional force on load y", "spring force", ...
    "unlatching motor force"];

writematrix(headers, csvFilePath);%creates headers for output file
writematrix(sol, csvFilePath, 'WriteMode', 'append');% addes actual data to csv file

end