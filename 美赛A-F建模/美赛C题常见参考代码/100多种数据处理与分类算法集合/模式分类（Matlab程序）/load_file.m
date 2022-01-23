function [patterns, targets, distribution_parameters] = load_file(filename)

%Load a file with either data points and/or distribution parameters

patterns = []; targets = []; distribution_parameters = [];

if isempty(findstr(filename,'.mat')),
    filename = [filename '.mat'];
end

if (~isempty(dir(filename)))
    load (filename)
    
    if (exist('features'))
        disp('Old data type. Changing the matrix called "features" into a matrix called "patterns". Please save the data once it is transformed');
        patterns = features;
    end
    
    if ((isempty(patterns) | isempty(targets))) & isempty(distribution_parameters) & ~exist('m0')
        error('No patterns, targets, or distribution parameters found in this file')
    end
    
    if (length(unique(targets)) ~= 2)
        error('The GUI requires two classes of data to operate. Use the text based interface for multiclass or single class data.')
    end
    
    hm = findobj('Tag', 'Messages');
    st= '';
    if (~isempty(patterns)),
        st = ['File loaded. Found ' num2str(size(targets,2)) ' data points'];
        
        if (size(patterns,1) > 2),
            %More than two dimensions in the data
            h = feature_selection;
            waitfor(h, 'UserData',1)
            
            h1		 = findobj(h, 'Tag', 'txtHiddenMethod');
            chosen = get(h1, 'String');
            h1		 = findobj(h, 'Tag', 'txtHiddenParams');
            params = get(h1, 'String');
            
            if (~isempty(str2num(params))),
                params = str2num(params);
            end
            
            [patterns, targets] = feval(chosen, patterns, targets, params);
            
            close(h)
        end
        
        param = max(max(abs(patterns)));
        plot_scatter(patterns, targets) 
        axis([-param,param,-param,param])
    end
    
    if ~isempty(distribution_parameters)
		% try (distribution_parameters.p0)
		%     %Old format
		%     disp('Old parametric structure found. Transforming into the new format. Please save the data once it is transformed');
		%     dp(1).mu    = distribution_parameters.m0;
		%     dp(2).mu    = distribution_parameters.m1;
		%     dp(1).sigma = distribution_parameters.s0;
		%     dp(2).sigma = distribution_parameters.s1;
		%     dp(1).w     = distribution_parameters.w0;
		%     dp(2).w     = distribution_parameters.w1;
		%     dp(1).p     = distribution_parameters.p0;
		%     dp(2).p     = 1 - distribution_parameters.p0;
		%     for i = 1:2,
		%         for j = 1:length(dp(i).w)
		%             dp(i).type(j,:) = cellstr('Gaussian');
		%         end
		%     end
		%     distribution_parameters = dp;
		% catch    
		% end
        
        n0 = size(distribution_parameters(1).sigma,1);
        n1 = size(distribution_parameters(2).sigma,1);
        if (~isempty(st)),
            st = [st ', '];
        else
            st = 'Found ';
        end
        st = [st num2str(n0) ' Gaussians for class 0 and ' num2str(n1) ' Gaussians for class 1.'];         
    end
    
    set(hm,'String',st)  
else
    hm = findobj('Tag', 'Messages');
    set(hm,'String','File not found.')
end 

if (nargout == 1),
    patterns = distribution_parameters;
end
