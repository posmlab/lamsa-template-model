function output_directory = create_output_directory()
curr_dir = pwd;
dateString = string(datetime);
cleanDateString = regexprep(dateString, " ", "_");
cleanDateString = regexprep(cleanDateString, ":", "_");
output_directory = fullfile('../output',cleanDateString);
if (~exist(fullfile(output_directory),'file'))
    mkdir(fullfile(output_directory));
end
cd(output_directory);
output_directory = pwd;
cd(curr_dir);