function [fh,ah,ch] = surf_custom(varargin)

p = parseInputs(varargin);

nAx = length(p.data1);         % How many axes?
if isempty(p.dim)              % Dimensions of subplots
    rows = round( sqrt(nAx) );
    cols = ceil(  sqrt(nAx) ); % More columns than rows by default
else
    rows = p.dim(1);
    cols = p.dim(2);
end

if     nAx == 1; fh = figure('Name',p.name,'Units','normalized','Position', [0, 0, 0.5, 0.5]); % Quarter screen
elseif nAx == 2; fh = figure('Name',p.name,'Units','normalized','Position', [0, 0, 1  , 0.5]); % Half screen
else;            fh = figure('Name',p.name,'Units','normalized','Position', [0, 0, 1  , 1  ]); % Full screen
end

% ah = tight_subplot(rows, cols, [0.08,0.08], [0.1,0.05], [0.04,0.05]);
ah = tight_subplot(rows, cols, [0.08,0.12], [0.1,0.12], [0.04,0.12]);
if strcmp(p.order,'rows')
    ah = reshape(ah,cols,rows);
    ah = ah.'; % Transpose the axes to switch the order
    ah = reshape(ah,rows*cols,1);
end

for ii = 1 : nAx
    
    axes(ah(ii));          % Make next axes the current axes
    pos = ah(ii).Position; % Save position information for later
    
    % Plot the surf image
    switch p.nData
        case 1; h = surf(p.data1{ii});
        case 2; h = surf(p.data1{ii}, p.data2{ii});
        case 3; h = surf(p.data1{ii}, p.data2{ii}, p.data3{ii});
        case 4; h = surf(p.data1{ii}, p.data2{ii}, p.data3{ii}, p.data4{ii});
    end
    
    % Choose colourmap
    colormap(p.color{ii});
    
    % Colourmap limits
    if isnumeric(p.clim{ii})
        caxis(p.clim{ii});
    end
    
    % Add colourbar
    if strcmp(p.colorbar{ii},'on')
        ch                = colorbar;
        ch.Label.String   = p.colorbarlabel{ii};
        ch.Label.FontSize = p.labelfontsize;
    end
    
    % Add labels with relevant font sizes
    set(gca,'FontSize',p.axisfontsize);
    title( p.title{ii} ,'FontSize',p.titlefontsize);
    xlabel(p.xlabel{ii},'FontSize',p.labelfontsize);
    ylabel(p.ylabel{ii},'FontSize',p.labelfontsize);
    
    % Misc
    set(h,'LineStyle',p.linestyle);
    view(p.view);
    axis tight;
    box on;
    
    % Change axes scales
    set(gca,'XScale',p.xscale{ii});
    set(gca,'YScale',p.yscale{ii});
    set(gca,'ZScale',p.zscale{ii});
    
    % XTick rules
    if strcmp(p.xtick{ii},'custom')
        xl = xlim;
        xticks(xl(1) : round((xl(2)-xl(1))/5,1,'significant') : xl(2))
    elseif isnumeric(p.xtick{ii})
        xticks(p.xtick{ii});
    end
    
    % YTick rules
    if strcmp(p.ytick{ii},'custom')
        yl = ylim;
        yticks(yl(1) : round((yl(2)-yl(1))/5,1,'significant') : yl(2))
    elseif isnumeric(p.ytick{ii})
        yticks(p.ytick{ii});
    end
    
    % ZTick rules
    if strcmp(p.ztick{ii},'custom')
        zl = zlim;
        zticks(zl(1) : round((zl(2)-zl(1))/5,1,'significant') : zl(2))
    elseif isnumeric(p.ztick{ii})
        zticks(p.ztick{ii});
    end
    
    % Lastly, adjust colourbar and label position
    if strcmp(p.colorbar{ii},'on')
        ch.Position(1) = pos(1)+pos(3)+0.002;
        ch.Label.Position(1) = 2.7;
    end
end

    function p = parseInputs(inputs)
        
        % Grab the inputs which should all be cells (passed in using {})
        p   =   inputParser;
        p.addParameter('data1'        , {0}         , @(x) iscell(x));
        p.addParameter('data2'        , {0}         , @(x) iscell(x));
        p.addParameter('data3'        , {0}         , @(x) iscell(x));
        p.addParameter('data4'        , {}         , @(x) iscell(x));
        p.addParameter('xlabel'       , {''}       , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('ylabel'       , {''}       , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('title'        , {''}       , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('name'         , ''         , @ischar);
        p.addParameter('color'        , {'default'}, @(x) iscell(x) && ischar(x{1}));
        p.addParameter('colorbar'     , {'on'}     , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('colorbarlabel', {''}       , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('xscale'       , {'linear'} , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('yscale'       , {'linear'} , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('zscale'       , {'linear'} , @(x) iscell(x) && ischar(x{1}));
        p.addParameter('xlim'         , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('xtick'        , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('ylim'         , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('ytick'        , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('zlim'         , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('ztick'        , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('clim'         , {'default'}, @(x) iscell(x) && (isnumeric(x{1}) || ischar(x{1})));
        p.addParameter('dim'          , []         , @(x) isnumeric(x) && length(x)==2);
        p.addParameter('axisfontsize' , 12         , @isscalar);
        p.addParameter('labelfontsize', 16         , @isscalar);
        p.addParameter('titlefontsize', 16         , @isscalar);
        p.addParameter('view'         , 2          , @isnumeric);
        p.addParameter('order'        , 'columns'  , @ischar   );
        p.addParameter('linestyle'    , 'none'     , @ischar);
        p.parse(inputs{:});
        p = p.Results;
        
        % How many surf inputs?
        if     isfield(p,'data4') && ~isempty(p.data4); p.nData = 4; 
        elseif isfield(p,'data3') && ~isempty(p.data3); p.nData = 3;
        elseif isfield(p,'data2') && ~isempty(p.data2); p.nData = 2;
        else;                                           p.nData = 1;
        end
    
        % If only one set of these inputs is supplied, copy it
        if length(p.data1) > 1 && length(p.xlabel) == 1
            for jj = 2:length(p.data1); p.xlabel{jj} = p.xlabel{1}; end
        end
        if length(p.data1) > 1 && length(p.ylabel) == 1
            for jj = 2:length(p.data1); p.ylabel{jj} = p.ylabel{1}; end
        end
        if length(p.data1) > 1 && length(p.title) == 1
            for jj = 2:length(p.data1); p.title{jj} = p.title{1}; end
        end
        if length(p.data1) > 1 && length(p.color) == 1
            for jj = 2:length(p.data1); p.color{jj} = p.color{1}; end
        end
        if length(p.data1) > 1 && length(p.colorbar) == 1
            for jj = 2:length(p.data1); p.colorbar{jj} = p.colorbar{1}; end
        end
        if length(p.data1) > 1 && length(p.colorbarlabel) == 1
            for jj = 2:length(p.data1); p.colorbarlabel{jj} = p.colorbarlabel{1}; end
        end
        if length(p.data1) > 1 && length(p.xscale) == 1
            for jj = 2:length(p.data1); p.xscale{jj} = p.xscale{1}; end
        end
        if length(p.data1) > 1 && length(p.yscale) == 1
            for jj = 2:length(p.data1); p.yscale{jj} = p.yscale{1}; end
        end
        if length(p.data1) > 1 && length(p.zscale) == 1
            for jj = 2:length(p.data1); p.zscale{jj} = p.zscale{1}; end
        end
        if length(p.data1) > 1 && length(p.xlim) == 1
            for jj = 2:length(p.data1); p.xlim{jj} = p.xlim{1}; end
        end
        if length(p.data1) > 1 && length(p.ylim) == 1
            for jj = 2:length(p.data1); p.ylim{jj} = p.ylim{1}; end
        end
        if length(p.data1) > 1 && length(p.zlim) == 1
            for jj = 2:length(p.data1); p.zlim{jj} = p.zlim{1}; end
        end
        if length(p.data1) > 1 && length(p.xtick) == 1
            for jj = 2:length(p.data1); p.xtick{jj} = p.xtick{1}; end
        end
        if length(p.data1) > 1 && length(p.ytick) == 1
            for jj = 2:length(p.data1); p.ytick{jj} = p.ytick{1}; end
        end
        if length(p.data1) > 1 && length(p.ztick) == 1
            for jj = 2:length(p.data1); p.ztick{jj} = p.ztick{1}; end
        end
        if length(p.data1) > 1 && length(p.clim) == 1
            for jj = 2:length(p.data1); p.clim{jj} = p.clim{1}; end
        end
        
    end



end