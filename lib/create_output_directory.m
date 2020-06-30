function output_directory = create_output_directory()
dateString = string(datetime);
cleanDateString = regexprep(dateString, " ", "_");
cleanDateString = regexprep(cleanDateString, ":", "_");
output_directory = fullfile('../output',cleanDateString);