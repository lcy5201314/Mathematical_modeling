function [patterns, targets, params, region] = click_points(region) 

%Manually enter points into the workspace
ax = region(1:4);

h   = findobj(findobj('Tag','classifier_GUI'),'Tag','txtNumberPointsPerClick');
if ~isempty(h),
    N   = str2num(get(h,'String'));
else
    N   = 20;
end
h   = findobj(findobj('Tag','classifier_GUI'),'Tag','Distribution');
if ~isempty(h),
    st  = get(h, 'String');
    val = get(h, 'Value');
    distribution = st(val);
else
    distribution = 'Gaussian';
end

hold on
axis(ax)
patterns    = [];
targets     = [];

distx			= (region(2) - region(1))/20;
disty			= (region(4) - region(3))/20;

%Print the text on figure:
t1 = text(region(1)+distx,region(3)+disty*3,'Class A  => Left  Mouse Click/Drag');
t2 = text(region(1)+distx,region(3)+disty*2,'Class B  => Right Mouse Click/Drag');
t3 = text(region(1)+distx,region(3)+disty,'To Stop  => Click any key on keyboard');

%Coordinates of figure for display:
h = findobj(gcbf,'Tag','classifier_GUI');
p = get(h,'Position');
wx1 = p(1);
wx2 = p(3);
wy1 = p(2);
wy2 = p(4);


d = 0.000001;

params(1).mu = [];
params(2).mu = [];

while 1
    %Take the x,y of patterns per one click 
    %k           = waitforbuttonpress;
    %point1      = get(h,'CurrentPoint');    % button down detected
    [x, y, button] = ginput(1);
    point1      = [x y];
    finalRect   = rbbox;                   % return figure units
    point2      = get(gca,'CurrentPoint');    % button up detected
    if (isempty(point1)),
        point1 = point2;
    end
    point1      = point1(1,1:2);              % extract x and y
    point2      = point2(1,1:2);
  
    mx   = (point1(1) + point2(1))/2;
    my   = (point1(2) + point2(2))/2;
    m    = [mx ; my];

    if (strcmp(distribution, 'Gaussian'))
        sx   = max(0.02, abs(point1(1)-point2(1))/5);
        sy   = max(0.02, abs(point1(2)-point2(2))/5);
        sigma = [sx 0 ; 0 sy]; 
        %Calculate the points:
        points = sigma * randn(2,N) + m * ones(1,N);
    else
        sx   = abs(point1(1)-point2(1));
        sy   = abs(point1(2)-point2(2));
        sigma = [sx, sy; 0, 0];
        %Calculate the points:
        points = (sigma(1,:)'*ones(1,N)) .* rand(2,N) + (m-sigma(1,:)'/2) * ones(1,N);
    end
    
    %Check the mouse output in order to set the class or exit :
    switch button
    case 1
        plot(points(1,:),points(2,:),'gx')
        patterns  = [patterns,points];
        targets   = [targets,ones(1,N)];
        in        = size(params(2).mu,1);
        params(2).mu(in+1,:)= m';
        params(2).sigma(in+1,:,:) = sigma;
        params(2).type(in+1,:) = cellstr(distribution);
    case 3
        plot(points(1,:),points(2,:),'bo')
        patterns  = [patterns,points];
        targets   = [targets,zeros(1,N)];
        in        = size(params(1).mu,1);
        params(1).mu(in+1,:)= m';
        params(1).sigma(in+1,:,:) = sigma;
        params(1).type(in+1,:) = cellstr(distribution);
    otherwise
        break
    end
   
    axis(ax)
end

in          = size(params(1).mu,1);
params(1).w = ones(in,1)/in;
in          = size(params(2).mu,1);
params(2).w = ones(in,1)/in;
params(1).p = length(params(1).w) / (length(params(1).w) + length(params(2).w));
params(2).p = length(params(2).w) / (length(params(1).w) + length(params(2).w));

region = [ax,100];

set(t1,'Visible','off')
set(t2,'Visible','off')
set(t3,'Visible','off')
hold off;
