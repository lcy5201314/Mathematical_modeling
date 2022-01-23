function plot_process(mu, plot_on)

%Plot the mu's during an algorithm's execution
%Inputs:
%   mu - The location of points
%   plot_on - 1-Plot centers, 2-Plot Voronoi centers, 3-Both
%
% No outputs

if (size(mu,1) > 2)
    %Can't plot more than two dimensions
    return;
end

h_mu        = findobj('UserData','mu');
h_voronoi   = findobj('UserData','process_Voronoi');
h_plot      = findobj(gca); h_plot=h_plot(1);

if (~isempty(h_mu)),
    delete(h_mu)
end
if (~isempty(h_voronoi)),
    delete(h_voronoi)
end

if (nargin == 1)
    return;
end

if (plot_on == 0)
    return;
end

if (~isempty(mu) & bitand(plot_on,1))
    h_mu = plot(mu(1,:),mu(2,:),'ro','LineWidth',2);
    set(h_mu,'UserData','mu');
end

if (~isempty(mu) & bitand(plot_on,2))
    colormap ('default')
    
    if (size(mu, 2) < 2)
        return;
    end
    
    %Compute the Voronoi regions
    N       = 100;
    Nmu     = size(mu,2);
    xlim    = get(h_plot, 'XLim');
    ylim    = get(h_plot, 'YLim');
    x       = linspace(xlim(1), xlim(2), N);
    y       = linspace(ylim(1), ylim(2), N);
    mx      = ones(N,1) * x;
    my      = y' * ones(1,N);
    flatxy  = [mx(:), my(:)]';

    dist    = zeros(Nmu, N^2);
    for i = 1:Nmu,
        dist(i,:) = (mu(1,i) - flatxy(1,:)).^2 + (mu(2,i) - flatxy(2,:)).^2;
    end
    [m, label]      = min(dist);
    [m, h_voronoi]  = contourf(x,y,reshape(label, N, N));
    for i = 1:size(h_voronoi,1)
        set(h_voronoi(i),'UserData','process_Voronoi');
    end
end

drawnow
    