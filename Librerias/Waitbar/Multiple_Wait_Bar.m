%% Provides a graphical series of wait bars in a single figure that doesn't interfere with updating multiple other figures
%:Inputs:
% - Progress (structure)
%       Field: Title - Name displayed above each waitbar
%       Field: Progress - Fraction (0-1) of progress
%       Field: Colour - Displayed colour of a waitbar
% - Figure Handle from Multiple_Wait_Bar Output (false if no figure able to be created)
%:Outputs:
% - Figure Handle
function Figure_Out = Multiple_Wait_Bar(Progress, Figure_Handle)
    %% Example Use:
    %AXIS TITLES
    %   Progress(1).Title = 'i';
    %   Progress(2).Title = 'j';
    %   Progress(3).Title = 'k';
    %BAR COLOURS
    %   Progress(1).Colour = 'r';
    %   Progress(2).Colour = 'g';
    %   Progress(3).Colour = 'b';
    %INITIATE FIGURE
    %   Progress_Figure = Multiple_Wait_Bar(Progress);
    %   for i = 0:0.2:1
    %       Progress(1).Progress = i;
    %       for j = 0:0.2:1
    %           Progress(2).Progress = j;
    %           for k = 0:0.1:1
    %               Progress(3).Progress = k;
    %DISPLAY PROGRESS
    %               Progress_Figure = Multiple_Wait_Bar(Progress, Progress_Figure);
    %           end
    %        end
    %   end
    %REMOVE VARIABLES FROM MEMORY
    %   clear('Progress_Figure','Progress');

    %% Validate Function Inputs
    if((nargin == 1) ||(nargin == 2))
        %validate input structure exists
        if ~isstruct(Progress)
            error('Expected Structure Input')
        else
            %get structure fieldnames
            Structure_Fieldnames = fieldnames(Progress);
            %check structure is empty
            if(isempty(Structure_Fieldnames))
                error('No Information contained within structure')
            else
                %get number of axis to plot
                NumAxes = length(Progress);
                %check if colour field exists
                if(~ isfield(Progress,'Colour'))
                    %enter empty title
                    [Progress(:).Colour] = deal('r');
                end
                %check if title field exists
                if(~ isfield(Progress,'Title'))
                    %enter empty title
                    [Progress(:).Title] = deal('');
                end
                %check progress field exists
                if(~ isfield(Progress,'Progress'))
                    %assume progress to be 0% completed
                    [Progress(:).Progress] = deal(0);
                end
                TitleAxes = {Progress(:).Title};
                StepAxes = [Progress.Progress];
            end
        end
    end
    %validate axis title data types
    if ~iscellstr(TitleAxes)
        error('Titles should be in a String')
    end
    %verify size of array
    SS = size(StepAxes,2);
    ST = size(TitleAxes,2);
    if ~isequal(SS,ST) && ~isequal(NumAxes,ST) && ~isequal(SS,NumAxes)
        error('Invalid Input Arguments');
    end
    %% check figure still exists
    Create_Figure = false;
    if(nargin >=2)
        %check figure exists
        if(~islogical(Figure_Handle))
            if(~ ishandle(Figure_Handle.figure))
                Create_Figure = true;
            end
        else
            Create_Figure = true;
            clear Figure_Handle;
        end
    end
    %% Create or Update Figure
    try
        %if figure requires creation
        if ((nargin < 2) || (Create_Figure))
            set(0, 'Units', 'pixel');
            screenSize = get(0,'ScreenSize');
            pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
            width = 360 * pointsPerPixel;
            height = NumAxes * 75 * pointsPerPixel;
            pos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];
            Figure_Handle.figure = figure('Position', pos, 'MenuBar', 'none',...
                'Numbertitle', 'off', 'Name', 'Please Wait..', 'Resize', 'off');
            axeswidth = 172;
            axesheight = 172/14;
            for i = 1:NumAxes
                R = max(0,min(100*Progress(i).Progress,100));
                X = 15;
                Y = 15+((NumAxes-i)*(pos(4)/(1.4*NumAxes)));
                axPos = [X Y axeswidth axesheight];
                Figure_Handle.Axeshandle(i).list = axes('Parent', Figure_Handle.figure, ...
                    'Box', 'on', 'Units', 'Points', 'Position', axPos);
                try
                    fill([0 R R 0], [0 0 1 1], Progress(i).Colour);
                catch
                    fill([0 R R 0], [0 0 1 1], 'r');
                end
                set(Figure_Handle.Axeshandle(i).list,'XLim',[0 100], 'YLim', [0 1], 'XTick', [], 'YTick', []);
                titlehand = get(Figure_Handle.Axeshandle(i).list, 'title');
                set(titlehand, 'string', cellstr(Progress(i).Title));
                drawnow;
            end
        %update figure    
        else
            for i = 1:NumAxes
                R = max(0,min(100*Progress(i).Progress,100));
                try
                    fill(Figure_Handle.Axeshandle(i).list, [0 R R 0], [0 0 1 1], Progress(i).Colour);
                catch
                    fill(Figure_Handle.Axeshandle(i).list, [0 R R 0], [0 0 1 1], 'r');
                end
                set(Figure_Handle.Axeshandle(i).list,'XLim',[0 100], 'YLim', [0 1], 'XTick', [], 'YTick', []);
                titlehand = get(Figure_Handle.Axeshandle(i).list,'title');
                set(titlehand,'string',cellstr(Progress(i).Title));
                drawnow;
                end
        end
    %% catch errors when drawing the waitbars (figure likely destroyed during process)
    catch
        %attempt to recreate the figure
        try
            set(0, 'Units', 'pixel');
            screenSize = get(0,'ScreenSize');
            pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
            width = 360 * pointsPerPixel;
            height = NumAxes * 75 * pointsPerPixel;
            pos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];
            Figure_Handle.figure = figure('Position', pos, 'MenuBar', 'none',...
                'Numbertitle', 'off', 'Name', 'Please Wait..', 'Resize', 'off');
            axeswidth = 172;
            axesheight = 172/14;
            for i = 1:NumAxes
                R = max(0,min(100*Progress(i).Progress,100));
                X = 15;
                Y = 15+((NumAxes-i)*(pos(4)/(1.4*NumAxes)));
                axPos = [X Y axeswidth axesheight];
                Figure_Handle.Axeshandle(i).list = axes('Parent', Figure_Handle.figure, ...
                    'Box', 'on', 'Units', 'Points', 'Position', axPos);
                try
                    fill([0 R R 0], [0 0 1 1], Progress(i).Colour);
                catch
                    fill([0 R R 0], [0 0 1 1], 'r');
                end
                set(Figure_Handle.Axeshandle(i).list,'XLim',[0 100], 'YLim', [0 1], 'XTick', [], 'YTick', []);
                titlehand = get(Figure_Handle.Axeshandle(i).list, 'title');
                set(titlehand, 'string', cellstr(Progress(i).Title));
                drawnow;
            end
        catch
            disp("Unable to create or update figure");
        end
    end
    %% Automatically close figure on completion of all bars
    if(mean([Progress.Progress]) == 1)
        if(ishandle(Figure_Handle.figure))
            close(Figure_Handle.figure);
        end
        clear Figure_Handle;
        Figure_Handle = false;
    end
    %% If an argument out is requested; provide figure handle to waitbar
    if nargout > 0
        Figure_Out = Figure_Handle;
    end
end