function show_algorithms (type, show_details)

% Specify possible classification algorithms and their details
%
% Inputs:
%   type            - Can be either classification, preprocessing, or feature_selection
%   show details    - (Optional) will show the necessary inputs and the default values for the algorithms

if (nargin == 1)
    show_details = 0;
end

switch lower(type)
case 'classification'
    fn = 'Classification.txt';
case 'preprocessing'
    fn = 'Preprocessing.txt';
case 'feature_selection'
    fn = 'Feature_selection.txt';
otherwise
    error('Unknown type.');
end

algorithms = read_algorithms(fn);

%Determine the width of the fields
Nname       = 10;
Ncaption    = 6;
for i = 1:length(algorithms)
    Nname       = max(Nname, length(algorithms(i).Name));
    Ncaption    = max(Ncaption, length(algorithms(i).Caption));
end
Nname    = Nname + 1;
Ncaption = Ncaption + 1;

%Show the fields
s = ['ALGORITHM', blanks(Nname - 9), 'INPUTS', blanks(Ncaption - 6), 'DEFAULT'];
disp(s)
disp(char(ones(1,length(s))*45))
for i = 1:length(algorithms)
    s = algorithms(i).Name;
    if (show_details)
        %Add spaces as needed
        s = [s, blanks(Nname - length(algorithms(i).Name))];
        
        if (strcmp(deblank(algorithms(i).Caption), ''))
            s = [s, 'None', blanks(Ncaption - 4)];
        else
            s = [s, algorithms(i).Caption];
            %Add spaces as needed
            s = [s, blanks(Ncaption - length(algorithms(i).Caption))];
        end
        
        s = [s, algorithms(i).Default];
    end
    disp(s)
end
