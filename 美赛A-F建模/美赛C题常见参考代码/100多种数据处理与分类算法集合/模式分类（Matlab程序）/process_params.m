function [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10] = process_params(params)

%This function recieves a parameter vector and returns it's components
[p1, p2, p3, p4, p5, p6, p7, p8, p9, p10] = deal([]);

if isnumeric(params),
    %The parameter vector is completely numeric. This is much simpler
    for i = 1:length(params),
        eval(['p' num2str(i) ' = ' num2str(params(i)) ';']);
    end
    Nassigned = length(params);   
else
    comma_loc       = findstr(params,',');
    Nassigned		 = 0;
    
    %First, locate any sub-parameter vectors
    bracket_loc1	= findstr(params,'[');
    bracket_loc2	= findstr(params,']');
    if (length(bracket_loc1) ~= length(bracket_loc2)),
        error('The number of opening and closing brackets in the parameter vectors is not equal!')
    end
    
    switch length(bracket_loc1),
    case 0
        error('No opening brackets in the parameter vector!')
    case 1
        if ((bracket_loc1 == 1) & (bracket_loc2 == length(params))),
            %This are just the usual opening and closing brackets
            params			= [',', params(2:end-1), ','];
            comma_loc      = findstr(params,',');
        else
            error('There are erroneous square brackets in the parameter vector!')
        end
    otherwise
        if ((bracket_loc1(1) == 1) & (bracket_loc2(end) == length(params))),      
            params			= [',', params(2:end-1), ','];
            bracket_loc1	= bracket_loc1(2:end);
            bracket_loc2	= bracket_loc2(1:end-1);
            comma_loc      = findstr(params,',');
        end
        for i = 1:length(bracket_loc1),
            %Find the next comma after the closing bracket
            close_loc	= bracket_loc2(i);
            open_loc		= bracket_loc1(i);
            close_comma = find(comma_loc < open_loc);
            close_comma = close_comma(end);
            val         = params(open_loc:close_loc);
            eval(['p' num2str(close_comma) ' = val;']);
            Nassigned	= Nassigned + 1;
        end
        
        %Now remove these sub-parameters
        for i = 1:length(bracket_loc1),
            bracket_loc1	= findstr(params,'['); bracket_loc1 = bracket_loc1(1);
            bracket_loc2	= findstr(params,']'); bracket_loc2 = bracket_loc2(1);
            if (bracket_loc1 ~= bracket_loc2 - 1),
                params			= params([1:bracket_loc1-1,bracket_loc2+1:end]);
            end
        end
        
    end
    
    %Now assign the remaining parameters to p's
    comma_loc      = findstr(params,',');
    for i = 1:length(comma_loc)-1,
        if (eval(['isempty(p' num2str(i) ')'])),
            tempp			= params(comma_loc(i)+1:comma_loc(i+1)-1);
            eval(['p' num2str(i) ' = ' deblank(tempp) ';']);
            Nassigned	= Nassigned + 1;
        end
    end
    
end

if (Nassigned ~= nargout),
    warning('Not all parameters in the parameter vector were read by the algorithm')
end
