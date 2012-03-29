function [d,xi,yi] = roipoly(x,y,a,xi,yi)
%ROIPOLY Define polygonal region of interest (ROI).
%	[BW,Xi,Yi] = ROIPOLY prompts you to define a region of interest 
%	using the mouse.  The region of interest is returned in the 
%	binary image BW which has the same size as the current image.  
%	BW contains 0's outside the region of interest and 1's inside.  
%	Also returned are the vertices of the polygon in Xi and Yi.
%	The ROI is defined by clicking in the image to define a closed
%	polygon. Press RETURN or right mouse button (shift-click on
%	the Macintosh) when done.
%
%	BW = ROIPOLY(A,Xi,Yi) returns a binary image BW that contains 1's 
%	inside the polygon defined by the vertices (x,y) and 0's
%	outside.
%
%	BW = ROIPOLY(X,Y,A,Xi,Yi) accounts for non-default axis limits.
%
%	Note:  the algorithm used by this function uses region
%	boundaries that lie entirely along edges between pixels.
%	For example, to select the rectangular region containing the
%	pixels A(a:b, c:d), set XI = [a-.5 a-.5 c+.5 c+.5] and
%	YI = [b-.5 d+.5 d+.5 b-.5] (assuming default axis limits).
%
%	See also ROICOLOR.

%	Reference:  For a discussion on region filling using boundaries
%	that lie entirely between pixels, see the technical report
%	"A Region Fill Algorithm Based on Crack Code Boundaries," by
%	Robert Grzeszczuk, Steven L. Eddins, and Thomas DeFanti,
%	Department of Electrical Engineering and Computer Science,
%	University of Illinois at Chicago, Chicago, Illinois 60680-4348,
%	UIC-EECS-92-6, 1992.

%	Clay M. Thompson 1-28-93
%	Revised Steven L. Eddins 14 March 1994
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.23 $  $Date: 1994/12/16 17:59:12 $


if nargin==0, % Get information from the current figure
  [x,y,a,hasimage] = getimage;
  if ~hasimage,
    error('The current figure must contain an image to use ROIPOLY.');
  end
%  [xi,yi] = getline(gcf); % Get rect info from the user.
  [xi,yi] = mpi_ginput;    % or no button presses.

elseif nargin==1 | nargin==2 | nargin==4,
  error('Wrong number of input arguments.');

elseif nargin==3,
  yi = a(:); a = x; xi = y(:); 
  [num_rows,num_cols] = size(a);
  x = [1 num_cols]; y = [1 num_rows];
else
  xi = xi(:); yi = yi(:);
  x = [x(1) x(prod(size(x)))];
  y = [y(1) y(prod(size(y)))];
end

[num_rows, num_cols] = size(a);

if length(xi)~=length(yi)
  error('XI and YI must be the same length.'); 
end
if length(xi)==0
  error('XI and YI can''t be empty.'); 
end

% Make sure polygon is closed.
xi = [xi;xi(1)]; yi = [yi;yi(1)];
% Transform xi,yi into pixel coordinates.
roix = axes2pix(num_cols, x, xi);
roiy = axes2pix(num_rows, y, yi);

% Round input vertices to the nearest 0.5 and then add 0.5.
roix = floor(roix + 1);
roiy = floor(roiy + 1);
if ((min(roiy(:)) < 1) | (max(roiy(:)) > (num_rows + 1)))
  error('Y vector contains out-of-range elements.');
end
if ((min(roix(:)) < 1) | (max(roix(:)) > (num_cols + 1)))
  error('X vector contains out-of-range elements.');
end

% Test for special case:  rectangular ROI.
if (isrect(roix, roiy))
  xmin = min(roix);
  xmax = max(roix);
  ymin = min(roiy);
  ymax = max(roiy);
  d = zeros(num_rows, num_cols);
  d(ymin:(ymax-1), xmin:(xmax-1)) = ones(ymax-ymin, xmax-xmin);
  return;
end

% Initialize output matrix.  We need one extra row to begin with.
d = zeros(num_rows+1,num_cols);

num_segments = max(size(roix)) - 1;

% Process each segment.
for counter = 1:num_segments
  x1 = roix(counter);
  x2 = roix(counter+1);
  y1 = roiy(counter);
  y2 = roiy(counter+1);

  % We only have to do something with this segment if it is not vertical
  % or a single point.
  if (x1 ~= x2)

    % Compute an approximation to the segment drawn on an integer
    % grid.  Mark appropriate changes in the x direction in the
    % output image.
    [x,y] = intline(x1,x2,y1,y2);
    diffx = diff(x);
    dx_indices = find(diffx);
    if (x2 > x1)
      mark_val = 1;
    else
      mark_val = -1;
      dx_indices = dx_indices + 1;
    end
    d_indices = [y(dx_indices) (x(dx_indices)-1)] * [1; (num_rows+1)];
    d(d_indices) = d(d_indices) + mark_val(ones(size(d_indices)),1);
  end
    
end

% Now a cumulative sum down the columns will fill the region with 
% either 1's or -1's.
d = abs(cumsum(d));

% Get rid of that extra row and we're done!
d = d(1:num_rows,:);
function pixelx = axes2pix(dim, x, axesx)
%AXES2PIX Convert axes coordinates to pixel coordinates
%	PIXELX = AXES2PIX(DIM, X, AXESX) converts axes coordinates
%	(as returned by get(gca, 'CurrentPoint'), for example) into
%	pixel coordinates.  X should be the vector returned by
%	X = get(image_handle, 'XData') (or 'YData').  DIM is the
%	number of image columns for the x coordinate, or the number
%	of image rows for the y coordinate.

%	See also IMAGE.

%	Steven L. Eddins 14 March 1994
%	Copyright (c) 1994 by The MathWorks, Inc.
%	$Revision: 1.3 $  $Date: 1994/11/29 21:27:50 $

% Error checking on input arguments.
error_str = nargchk(3, 3, nargin);
if (~ isempty(error_str));
  error('There must be 3 input arguments.');
end
if (max(size(dim)) ~= 1)
  error('First argument must be a scalar.');
end
if (min(size(x)) > 1)
  error('X must be a vector.');
end

xfirst = x(1);
xlast = x(max(size(x)));

if (dim == 1)
  pixelx = axesx - xfirst + 1;
  return;
end
xslope = (dim - 1) / (xlast - xfirst);
if ((xslope == 1) & (xfirst == 1))
  pixelx = axesx;
else
  pixelx = xslope * (axesx - xfirst) + 1;
end
function varargout = getimage(varargin)
%GETIMAGE Get image data from axes.
%   A = GETIMAGE(H) returns the first image data contained in
%   the Handle Graphics object H.  H can be a figure, axes,
%   image, or texture-mapped surface.  A is identical to the
%   image CData; it contains the same values and is of the same
%   class (uint8 or double) as the image CData. If H is not an
%   image or does not contain an image or texture-mapped surface,
%   A is empty.
%
%   [X,Y,A] = GETIMAGE(H) returns the image XData in X and the
%   YData in Y. XData and YData are two-element vectors that
%   indicate the range of the x-axis and y-axis.
%
%   [...,A,FLAG] = GETIMAGE(H) returns an integer flag that
%   indicates the type of image H contains. FLAG is one of these
%   values:
%   
%       0   not an image; A is returned as an empty matrix
%
%       1   indexed image
%
%       2   intensity image with values in standard range ([0,1]
%           for double arrays, [0,255] for uint8 arrays,
%           [0,65535] for uint16 arrays)
%
%       3   intensity data, but not in standard range
%
%       4   RGB image
%
%   [...] = GETIMAGE returns information for the current axes. It
%   is equivalent to [...] = GETIMAGE(GCA).
%
%   Class Support
%   -------------
%   The output array A is of the same class as the image
%   CData. All other inputs and outputs are of class double.
%
%   Example
%   -------
%   After using imshow to display an image directly from a file,
%   use GETIMAGE to get the image data into the workspace.
%
%       imshow rice.tif
%       I = getimage;

%   Copyright 1993-2001 The MathWorks, Inc.  
%   $Revision: 5.18 $  $Date: 2001/01/18 15:29:07 $

% Subfunctions:
% - FINDIM
% - GET_IMAGE_INFO

him = findim(varargin{:});

[x,y,A,state] = get_image_info(him);

switch nargout
case 0
    % GETIMAGE(...)
    varargout{1} = A;
    
case 1
    % A = GETIMAGE(...)
    varargout{1} = A;
    
case 2
    % [A,FLAG] = GETIMAGE(...)
    varargout{1} = A;
    varargout{2} = state;
    
case 3
    % [x,y,A] = GETIMAGE(...)
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = A;
    
case 4
    % [x,y,A,FLAG] = GETIMAGE(...)
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = A;
    varargout{4} = state;
    
otherwise
    error('Too many output arguments.');
    
end

%----------------------------------------------------------------------
% Local Function: FINDIM
%----------------------------------------------------------------------
function him = findim(varargin)
%FINDIM Find image object.
%   HIM = FINDIM(H) searches for a valid Handle Graphics Image
%   object starting from the handle H and returns its handle in
%   HIM. H may be the handle of a Figure, Axes, or Image object.
%
%   If H is an Image object, FINDIM returns it.
%
%   If H is an Axes object, FINDIM searches H for Image objects.
%   If more than one Image object is found in the Axes, FINDIM
%   looks to see if one of the Images is the current object. If
%   so, FINDIM returns that Image. Otherwise, FINDIM returns the
%   highest Image in the stacking order.
%
%   If H is a Figure object, FINDIM searches H's current Axes.
%
%   HIM = FINDIM searchs the current Figure.

him = [];
if (nargin == 0)
    rootKids = get(0,'Children');
    if (~isempty(rootKids))
        figHandle = get(0,'CurrentFigure');
        figAxes = findobj(get(figHandle, 'Children'), 'flat', 'Type', 'axes');
        if (~isempty(figAxes))
            axHandle = get(figHandle, 'CurrentAxes');
            him = findim_in_axes(axHandle);
        end
    end
else
    % User specified a handle.
    h = varargin{1};
    h = h(1);
    if (~ishandle(h))
        error('Invalid handle H.');
    end
    switch get(varargin{1},'Type')
    case 'figure'
        figHandle = varargin{1};
        figAxes = findobj(get(figHandle, 'Children'), 'flat', 'Type', 'axes');
        if (~isempty(figAxes))
            axHandle = get(figHandle, 'CurrentAxes');
            him = findim_in_axes(axHandle);
        end
        
    case 'axes'
        axHandle = varargin{1};
        him = findim_in_axes(axHandle);
        
    case 'image'
        him = h;
        
    otherwise
        error('Input handle H must be a figure, axes, or image.');
        
    end
    
end

%----------------------------------------------------------------------
% Local Function: FINDIM_IN_AXES
%----------------------------------------------------------------------
function him = findim_in_axes(axHandle)

figHandle = get(axHandle, 'Parent');
% If the current object is a texture-mapped surface, use that.
currentObj = get(figHandle, 'CurrentObject');
if (~isempty(currentObj) & strcmp(get(currentObj,'type'),'surface') & ...
            strcmp(get(currentObj,'FaceColor'),'texturemap'))
    him = currentObj;
else
    him = findobj(axHandle, 'Type', 'image');
    if (length(him) > 1)
        % Found more than one image in the axes.
        % If one of the images is the current object, use it.
        % Otherwise, use the first image in the stacking order.
        if (isempty(currentObj))
            % No current object; use the one on top.
            him = him(1);
        else
            % If the current object is one of the images
            % we found, use it.
            idx = find(him == currentObj);
            if (isempty(idx))
                him = him(1);
            else
                him = him(idx);
            end
        end
    end
end
if (isempty(him))
    % Didn't find an image.  Is there a texturemapped surface we can use?
    him = findobj(axHandle, 'Type', 'surface', ...
            'FaceColor', 'texturemap');
    if (~isempty(him))
        him = him(1);
    end
end
            
            
%----------------------------------------------------------------------
% Local Function: GET_IMAGE_INFO
%----------------------------------------------------------------------
function [x,y,A,state] = get_image_info(him)

if (isempty(him))
    % We didn't find an image.
    x = [];
    y = [];
    A = [];
    state = 0;
    
elseif (strcmp(get(him, 'Type'), 'surface'))
    % We found a texturemapped surface object.
    A = get(him, 'CData');
    x = get(him, 'XData');
    y = get(him, 'YData');
    state = 2;
    
else
    % We did find an image.  Find out about it.
    userdata = get(him, 'UserData');
    cdatamapping = get(him, 'CDataMapping');
    x = get(him, 'XData');
    y = get(him, 'YData');
    A = get(him, 'CData');
    
    if ((ndims(A) == 3) & (size(A,3) == 3))
        % We have an RGB image
        state = 4;
        
    else
        % Not an RGB image
    
        if (isequal(cdatamapping,'direct'))
            % Do we have an indexed image or an old-style intensity
            % or scaled image?
            
            if (isequal(size(userdata), [1 2]))
                % We have an old-style intensity or scaled image.
                
                % How long is the colormap?
                N = size(get(get(get(him,'Parent'),'Parent'),'Colormap'),1);
                
                if (isequal(userdata, [0 1]))
                    % We have an old-style intensity image.
                    A = (A-1)/(N-1);
                    state = 2;
                    
                else
                    % We have an old-style scaled image.
                    A = (A-1)*((userdata(2)-userdata(1))/(N-1))+userdata(1);
                    state = 3;
                    
                end
                
            else
                % We have an indexed image.
                state = 1;
                
            end
            
        else
            % CDataMapping is 'scaled'
            
            hax = get(him, 'Parent');
            clim = get(hax, 'CLim');
            if ((isa(A,'double') & isequal(clim,[0 1])) | ...
                  (isa(A,'uint8') & isequal(clim,[0 255])) | ...
                  (isa(A,'uint16') & isequal(clim,[0 65535])))
                % We have an intensity image.
                state = 2;
                
            else
                % We have a scaled image.
                state = 3;
                
            end
        end
        
    end
    
end
function [x,y] = intline(x1, x2, y1, y2)
%INTLINE Integer-coordinate line drawing algorithm.
%	[X, Y] = INTLINE(X1, X2, Y1, Y2) computes an
%	approximation to the line segment joining (X1, Y1) and
%	(X2, Y2) with integer coordinates.  X1, X2, Y1, and Y2
%	should be integers.  INTLINE is reversible; that is,
%	INTLINE(X1, X2, Y1, Y2) produces the same results as
%	FLIPUD(INTLINE(X2, X1, Y2, Y1)).

%	Steven L. Eddins, October 1994
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.4 $  $Date: 1994/10/28 17:28:07 $

dx = abs(x2 - x1);
dy = abs(y2 - y1);

% Check for degenerate case.
if ((dx == 0) & (dy == 0))
  x = x1;
  y = y1;
  return;
end

flip = 0;
if (dx >= dy)
  if (x1 > x2)
    % Always "draw" from left to right.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (y2 - y1)/(x2 - x1);
  x = (x1:x2).';
  y = round(y1 + m*(x - x1));
else
  if (y1 > y2)
    % Always "draw" from bottom to top.
    t = x1; x1 = x2; x2 = t;
    t = y1; y1 = y2; y2 = t;
    flip = 1;
  end
  m = (x2 - x1)/(y2 - y1);
  y = (y1:y2).';
  x = round(x1 + m*(y - y1));
end
  
if (flip)
  x = flipud(x);
  y = flipud(y);
end
function r = isrect(xi, yi, arg3)
%ISRECT Determine if polygon is a rectangle.
%	ISRECT(XI, YI) returns 1 if (XI,YI) are the vertices of a
%	polygon that is a rectangle.  Otherwise it returns 0.
%
%	ISRECT(XI, YI, 'parallel') returns 1 if (XI, YI) are the
%	vertices of a polygon that is a rectangle oriented
%	parallel to the x and y axes.  Otherwise it returns 0.

%	Steven L. Eddins, 22 March 1994
%	Copyright (c) 1994 by The MathWorks, Inc.
%	$Revision: 1.5 $  $Date: 1995/01/11 20:43:59 $

error(nargchk(2,3,nargin));

xi = xi(:)';
yi = yi(:)';
N = size(xi,2);
if (N < 2)
  r = 1;
  return;
end
if ((xi(N) ~= xi(1)) | (yi(N) ~= yi(1)))
  % Close the polygon and add a point.
  xi = [xi  xi(1:2)];
  yi = [yi  yi(1:2)];
else
  xi = [xi xi(2)];
  yi = [yi yi(2)];
end

tol1 = max(max(abs(xi)), max(abs(yi))) * 10000 * eps;
tol2 = 10000*eps;

v = [diff(xi) ; diff(yi)];

% Remove duplicates
v(:,(abs(v(1,:)) < tol1) & (abs(v(2,:)) < tol1)) = [];
N = size(v,2);
if (N == 0)
  % Degenerate point rectangle.
  r = 1;
  return;
end

% Normalize the vectors
mags = sqrt(sum(v.^2));
v = v ./ [mags ; mags];

if (nargin == 3)
  if (any(prod(v) > tol2))
    % Check for non-horizontal or non-vertical vectors.
    r = 0;
    return;
  end
end

% Find dot products.
absdots = abs(sum(v(:,1:N-1) .* v(:,2:N)));

% Find right angles.
rights = absdots < tol2;
numRights = sum(rights);

% Find parallels.
parallels = abs(1-absdots) < tol2;

if (all(rights | parallels) & (rem(numRights,4) == 0))
  r = 1;
else
  r = 0;
end

