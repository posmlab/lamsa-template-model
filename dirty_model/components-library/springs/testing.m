
%path = '\Subclasses'

%files = dir(strcat(path,'*Spring.m'))

%files = dir('*Spring.m');

% files = dir(fullfile('Subclasses', '*.m'));
% 
% for file = files'
%     
%     filename = file.name;
%     [pathstr, name, ext] = fileparts(filename);
%     name = name
% end
% 
% addpath 'Subclasses';
% 
% 
% for property = properties('LinearSpring')'
%     
%     property = property(1);
%     
% end


[filepath, ~, ~] = fileparts(mfilename('fullpath'));
%path = filepath
%currentfolder = pwd
%files = dir('*.m');
%file = files(1)

cd(filepath)
cd ../;
cd ../;
cd ../;

newFile = pwd

%cd components-library;

%str = "yo"
% 
% files = dir();
% file = files(1)