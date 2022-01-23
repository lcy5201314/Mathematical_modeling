function FindParametersFunctions(action)

Nlines  = 6;
short_action = action(1:4);
param			 = str2num(action(end));
if ~isempty(param),
   hNone 	= findobj('Tag',['btnNone' num2str(param)]);
   hString	= findobj('Tag',['btnString' num2str(param)]);
   hNumber	= findobj('Tag',['btnNumeric' num2str(param)]);
end

switch short_action
case 'Star'
   %Start the parameter tuning
   
   %Find who to tune
   hDiv			= findobj('Tag','txtDivisions');
   divisions	= str2num(get(hDiv,'String'));
   hFolds		= findobj('Tag','txtFolds');
   Nfolds		= str2num(get(hFolds,'String'));
   tune			= ones(1,Nlines);
   
   for i = 1:Nlines,
	   hNumber	= findobj('Tag',['btnNumeric' num2str(i)]);
	   hString	= findobj('Tag',['btnString' num2str(i)]);
      if (get(hNumber, 'Value') == 1),
         tune(i) = divisions;
         hSpacing	= findobj('Tag',['popSpacing' num2str(i)]);
         hFrom		= findobj('Tag',['txtLowerLimit' num2str(i)]);
         hTo		= findobj('Tag',['txtUpperLimit' num2str(i)]);
         
         if (get(hSpacing,'Value') == 1),
            parameters(i).grid = linspace(str2num(get(hFrom,'String')),str2num(get(hTo,'String')),divisions);
         else
	         parameters(i).grid = logspace(log10(str2num(get(hFrom,'String'))),log10(str2num(get(hTo,'String'))),divisions);
         end
         
         parameters(i).string = '';
      else
         if (get(hString, 'Value') == 1),
            hString = findobj('Tag', ['txtString' num2str(i)]);
            parameters(i).string = get(hString,'String');
            parameters(i).grid = inf;
         else
            parameters(i).grid = inf;
            parameters(i).string = '';
         end
      end
   end
   
   numeric		= find(tune == divisions);
   Nnumeric		= length(numeric);
   best			= zeros(Nfolds, Nnumeric);
   
   %Prepare the data
   patterns		= evalin('base','patterns');
   targets		= evalin('base','targets');
   indices		= randperm(floor(size(patterns,2)/Nfolds)*Nfolds);
   indices		= reshape(indices,Nfolds,floor(size(patterns,2)/Nfolds));
   
   %Find the classifier
	h		      = findobj('Tag', 'Algorithm'); 
	val			= get(h,'Value');     	
	algorithm	= get(h,'String'); 	
	algorithm	= deblank(char(algorithm(val,:)));
   
   for f	= 1:Nfolds,
      %For each of the parameter options, find the classification results
      min_error = 1;
      disp(['Working on fold number ' num2str(f)])
      
      for i1 = 1:tune(1),
         for i2 = 1:tune(2),
            for i3 = 1:tune(3),
               for i4 = 1:tune(4),
                   for i5 = 1:tune(5),
                       for i6 = 1:tune(6),
                           %Build a parameter vector
                           Nvec	= zeros(1,Nnumeric);
                           param	= '[';
                           for j = 1:Nlines,
                               if isempty(parameters(j).string),
                                   if (length(parameters(j).grid) > 1),
                                       Nvec(find(numeric==j)) = eval(['parameters(j).grid(i' num2str(j) ')']);
                                       param	= [param, num2str(Nvec(find(numeric==j))) ','];
                                   end
                               else
                                   if isempty(str2num(parameters(j).string)),
                                       param = [param, '''', parameters(j).string, ''',' ];
                                   else
                                       param = [param, parameters(j).string, ','];
                                   end
                               end
                           end
                           
                           param = [param(1:end-1), ']'];
                           if ~isempty(str2num(param)),
                               param = str2num(param);
                           end
                           
                           %Build a classifier
                           train_in	= indices([1:f-1,f+1:Nfolds],:);
                           train_in = train_in(:);
                           test_in  = indices(f,:);
                           
                           test_targets	= feval(algorithm, patterns(:, train_in), targets(train_in), patterns(:, test_in), param);
                           test_err     = mean(test_targets ~= targets(test_in));                             
                           
                           if (test_err < min_error),
                               best(f,:) = Nvec;
                               min_error = test_err;
                           end
                       end
                   end
               end
            end	
         end
      end
   end
   
   %Finished on all the dimensions, so now find the optimal configuration
   best_params	= median(best')';
   
   %Form the best parameter vector
   param	= '[';
   for j = 1:Nlines,
      if isempty(parameters(j).string),
         if (length(parameters(j).grid) > 1),
            param	= [param, num2str(best_params(find(numeric==j))) ','];
         end
      else
         param = [param, '''', parameters(j).string ,''','];
      end
   end
   
   param = [param(1:end-1), ']'];
   
   %Put it inside the GUI
   uiwait(helpdlg(['The best training parameters found are: ', param], 'Parameter selection'))
   if ~isempty(str2num(param)),
      param = str2num(param);
   end
   
   if strcmp(get(findobj('Tag', 'txtAlgorithmParameters'),'Visible'),'on'),
	   set(findobj('Tag', 'txtAlgorithmParameters'),'String',param)
	else
	   set(findobj('Tag', 'txtAlgorithmParametersLong'),'String',param)
	end
   
   close
   
case 'None'     
   set(hNone, 'Value', 1);
   set(hString, 'Value', 0);
   set(hNumber, 'Value', 0);   
case 'Stri'
   set(hNone, 'Value', 0);
   set(hString, 'Value', 1);
   set(hNumber, 'Value', 0);
case 'Numb'
   set(hNone, 'Value', 0);
   set(hString, 'Value', 0);
   set(hNumber, 'Value', 1);
otherwise
   error('Unknown command')
end

