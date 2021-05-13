
%path = '\Subclasses'

%files = dir(strcat(path,'*Spring.m'))

%files = dir('*Spring.m');

files = dir(fullfile('Subclasses', '*.m'));

for file = files'
    
    filename = file.name;
    [pathstr, name, ext] = fileparts(filename);
    name = name
end

addpath 'Subclasses';


for property = properties('LinearSpring')'
    
    property = property(1);
    
end
