function algorithms = read_algorithms(FileName)

%This functions reads an algorithms file
%The file is a test file, where each line is as follows:
%<Algorithm name>;<Parameters name>;'<Default parameters>';<S/L>
%Where: S/L is wheather to display a short or a long parameter box. This field is optional
%The function returns a structure array built accordingly

fid = fopen(FileName,'r');

if (fid < 0),
    error(['The file ' FileName ' was not found. Cannot continue'])
end

all_file = fread(fid);
fclose(fid);

line_ends = find(all_file == 13);
line_ends = [-1 ; line_ends];
for i = 1:length(line_ends)-1,
    line = all_file(line_ends(i)+2:line_ends(i+1)-1)';
    delimiter = find(line == 64); %The @ sign
    if (length(delimiter) > 0),
        algorithms(i).Name = char(line(1:delimiter(1)-1));
        algorithms(i).Caption = char(line(delimiter(1)+1:delimiter(2)-1));
        if (length(delimiter) > 2),
            algorithms(i).Default = char(line(delimiter(2)+1:delimiter(3)-1));
            algorithms(i).Field   = char(line(delimiter(3)+1));
        else
            algorithms(i).Default = char(line(delimiter(2)+1:length(line)));
        end
    end
end
