function [img,H] = iradonKC_Single(varargin)
%iradonKC_Single compute image from k-space data acquired by
% radial readouts using filtered backprojection algorithm.
% Assumption:
% The center of k-space is sampled and corresponds to Nx/2+1 sample.
% 
% iradonKC_Single.m is based on iradon.m and has the same input/output
% interface. Exception: input data R is 2D or 3D array containing k-space
% data from single or multiple receiver coils (Nr-by-Np-by-Nc),
% where Nr - the number points in readout
%       Np - the munber of readouts (projections)
%       Nc - the number of coils
%
%            The description of iradon.m
%===============================================================
%   I = iradon(R,THETA) reconstructs the image I from projection 
%   data in the 2-D array R.  The columns of R are parallel beam 
%   projection data.  IRADON assumes that the center of rotation
%   is the center point of the projections, which is defined as
%   ceil(size(R,1)/2).
%
%   THETA describes the angles (in degrees) at which the projections    
%   were taken.  It can be either a vector containing the angles or 
%   a scalar specifying D_theta, the incremental angle between 
%   projections. If THETA is a vector, it must contain angles with 
%   equal spacing between them.  If THETA is a scalar specifying   
%   D_theta, the projections are taken at angles THETA = m * D_theta; 
%   m = 0,1,2,...,size(R,2)-1.  If the input is the empty matrix
%   ([]), D_theta defaults to 180/size(R,2).   
%
%   IRADON uses the filtered backprojection algorithm to perform
%   the inverse Radon transform.  The filter is designed directly 
%   in the frequency domain and then multiplied by the FFT of the 
%   projections.  The projections are zero-padded to a power of 2 
%   before filtering to prevent spatial domain aliasing and to 
%   speed up the FFT.
%
%   I = IRADON(R,THETA,INTERP,FILTER,D,N) specifies parameters 
%   to use in the inverse Radon transform.  You can specify
%   any combination of the last four arguments.  IRADON uses 
%   default values for any of these arguments that you omit.
%
%   INTERP specifies the type of interpolation to use in the   
%   backprojection.  The available options are listed in order
%   of increasing accuracy and computational complexity:
%
%      'nearest' - nearest neighbor interpolation 
%      'linear'  - linear interpolation (default)
%      'spline'  - spline interpolation
%
%   FILTER specifies the filter to use for frequency domain filtering.  
%   FILTER is a string that specifies any of the following standard 
%   filters:
% 
%   'Ram-Lak'     The cropped Ram-Lak or ramp filter (default).  The    
%                 frequency response of this filter is |f|.  Because 
%                 this filter is sensitive to noise in the projections, 
%                 one of the filters listed below may be preferable.   
%   'Shepp-Logan' The Shepp-Logan filter multiplies the Ram-Lak
%                 filter by a sinc function.
%   'Cosine'      The cosine filter multiplies the Ram-Lak filter 
%                 by a cosine function.
%   'Hamming'     The Hamming filter multiplies the Ram-Lak filter 
%                 by a Hamming window.
%   'Hann'        The Hann filter multiplies the Ram-Lak filter by 
%                 a Hann window.
%   
%   D is a scalar in the range (0,1] that modifies the filter by 
%   rescaling its frequency axis.  The default is 1.  If D is less 
%   than 1, the filter is compressed to fit into the frequency range  
%   [0,D], in normalized frequencies; all frequencies above D are set  
%   to 0.
% 
%   N is a scalar that specifies the number of rows and columns in the 
%   reconstructed image.  If N is not specified, the size is determined   
%   from the length of the projections:
%
%       N = 2*floor(size(R,1)/(2*sqrt(2)))
%
%   If you specify N, IRADON reconstructs a smaller or larger portion of 
%   the image, but does not change the scaling of the data.  
% 
%   If the projections were calculated with the RADON function, the 
%   reconstructed image may not be the same size as the original 
%   image.  
%
%   [I,H] = iradon(...) returns the frequency response of the filter
%   in the vector H.
%
%   Class Support
%   -------------
%   All input arguments must be of class double.  The output arguments are
%   of class double.
%
%   Example
%   -------
%       P = phantom(128);
%       R = radon(P,0:179);
%       I = iradon(R,0:179,'nearest','Hann');
%       imshow(P); figure; imshow(I);
%
%   See also RADON, PHANTOM.


[p,theta,filter,d,interp,N] = parse_inputs(varargin{:});

p = single(p);
theta = single(theta);


% Design the filter
len = size(p,1);   
H = designFilter(filter, len, d);
% EGK: The number of projections is finite. Therefore, non-zero weight for the central k-space point
H(1) = H(2)/8;

% EGK: No zero padding
%%p(length(H),1)=0;  % Zero pad projections 

% EGK: Projection data is k-space data for radial sampled MRI.
%%p = fft(p);    % p holds fft of projections

% EGK implemented vectorized version of filtering opeartion
p = fftshift(p,1);
Nc = size(p,3);

if Nc == 1
   try p = single(p.*repmat(H,1,size(p,2)));
   catch
     disp('Will crash off by one if input is odd. Need to have projections even length!')
     whos H
     whos p
     return
   end   %evrd 12/06
else
   try p = single(p.*repmat(H,[1,size(p,2),Nc]));
   catch
     disp('Will crash off by one if input is odd. Need to have projections even length!')
     return
   end
end

% EGK deleted real operator because for MRI we are working with complex data.
%%p = real(ifft(p));     % p is the filtered projections
p = fftshift(ifft(p),1);
clear i;

% EGK: No truncation, because no zero padding
%%p(len+1:end,:) = [];   % Truncate the filtered projections

if Nc == 1
   img = zeros(N,'single');        % Allocate memory for the image.
else
   img = zeros(N,N,Nc,'single');
end

% Define the x & y axes for the reconstructed image so that the origin
% (center) is in the spot which RADON would choose.
%%xax = (1:N)-ceil(N/2);
% EGK change the centre of projections to support Siemens implementation
 xax = (1:N)-(ceil(N/2)+1);
% xax = (1:N)-(N+1)/2;

x = single(repmat(xax, N, 1));    % x coordinates, the y coordinates are rot90(x)
y = single(rot90(x));

costheta = single(cos(theta));
sintheta = single(sin(theta));
%%ctrIdx = ceil(len/2);     % index of the center of the projections
% EGK change the center of the projections to support Siemens implementation
ctrIdx = ceil(len/2)+1;
%ctrIdx = (len+1)/2;

% Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
% quantity is greater than the length of the projections
imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.
if size(p,1) < imgDiag 
   rz = imgDiag - size(p,1);  % how many rows of zeros
   p = [zeros(ceil(rz/2),size(p,2)); p; zeros(floor(rz/2),size(p,2))];
   ctrIdx = ctrIdx+ceil(rz/2);
end

% Backprojection - vectorized in (x,y), looping over theta
if strcmp(interp, 'nearest neighbor')
   for m=1:length(theta)   
%%%      t = round(x*costheta(m) + y*sintheta(m));
      t = round(ctrIdx+x*costheta(m) + y*sintheta(m));
      for c=1:Nc
         proj = p(:,m,c);
%%%         img(:,:,c) = img(:,:,c) + proj(t+ctrIdx);
         img(:,:,c) = img(:,:,c) + proj(t);
      end
   end
elseif strcmp(interp, 'linear')
   for m=1:length(theta)
%%%      t = x.*costheta(m) + y.*sintheta(m);
      t = ctrIdx+x.*costheta(m) + y.*sintheta(m);
      a = floor(t);
      
      for c=1:Nc
         proj = p(:,m,c);
%%%         img(:,:,c) = img(:,:,c) + proj(a+ctrIdx) + (t-a).*(proj(a+1+ctrIdx)-proj(a+ctrIdx));
         img(:,:,c) = img(:,:,c) + proj(a) + (t-a).*(proj(a+1)-proj(a));
      end
   end
elseif strcmp(interp, 'spline')
   for m=1:length(theta)
%%%      taxis = (1:size(p,1)) - ctrIdx;
%%%      t = x.*costheta(m) + y.*sintheta(m);
      %taxis = (1:size(p,1));
      %t = ctrIdx + x.*costheta(m) + y.*sintheta(m);
      taxis = (1:size(p,1)) - floor(ctrIdx);
      t = ctrIdx - floor(ctrIdx) + x.*costheta(m) + y.*sintheta(m);
      for c=1:Nc
         proj = p(:,m,c);
         projContrib = interp1(taxis,proj,t(:),'*spline');
         img(:,:,c) = img(:,:,c) + reshape(projContrib,N,N);
      end
   end
end

img = img*pi/(2*length(theta));


%%%
%%%  Sub-Function:  designFilter
%%%

function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be 
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter 
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


% EGK: filter length is equal to projection length. No zero padding.
%%order = max(64,2^nextpow2(2*len));
order = len;

% First create a ramp filter - go up to the next highest
% power of 2.

filt = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch filter
case 'ram-lak'
   % Do nothing
case 'shepp-logan'
   % be careful not to divide by 0:
   filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
case 'cosine'
   filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
case 'hamming'  
   filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
case 'hann'
   filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
otherwise
   error('Invalid filter selected.');
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter


%%%
%%%  Sub-Function:  parse_inputs
%%%

function [p,theta,filter,d,interp,N] = parse_inputs(varargin);
%  Parse the input arguments and retun things
%
%  Inputs:   varargin -   Cell array containing all of the actual inputs
%
%  Outputs:  p        -   Projection data
%            theta    -   the angles at which the projections were taken
%            filter   -   string specifying filter or the actual filter
%            d        -   a scalar specifying normalized freq. at which to crop 
%                         the frequency response of the filter
%            interp   -   the type of interpolation to use
%            N        -   The size of the reconstructed image

if nargin<2
   error('Invalid input arguments.');
end

p = varargin{1};
theta = pi*varargin{2}/180;

% Default values
N = 0;                 % Size of the reconstructed image
d = 1;                 % Defaults to no cropping of filters frequency response
filter = 'ram-lak';    % The ramp filter is the default
interp = 'linear';     % default interpolation is linear
string_args = {'nearest neighbor', 'linear', 'spline', ...
      'ram-lak','shepp-logan','cosine','hamming', 'hann'};

for i=3:nargin
   arg = varargin{i};
   if ischar(arg)
      idx = strmatch(lower(arg),string_args);
      if isempty(idx)
         error(['Unknown input string: ' arg '.']);
      elseif prod(size(idx)) > 1
         error(['Ambiguous input string: ' arg '.']);
      elseif prod(size(idx)) == 1
         if idx <= 3   % It is the interpolatio
            interp = string_args{idx};
         elseif (idx > 3) & (idx <= 8)
            filter = string_args{idx};
         end
      end
   elseif prod(size(arg))==1
      if arg <=1
         d = arg;
      else
         N = arg;
      end
   else
      error('Invalid input parameters');
   end
end

% If the user didn't specify the size of the reconstruction, so 
% deduce it from the length of projections
if N==0    
   N = 2*floor( size(p,1)/(2*sqrt(2)) );  % This doesn't always jive with RADON
end

% for empty theta, choose an intelligent default delta-theta
if isempty(theta)
   theta = pi / size(p,2);
end

% If the user passed in delta-theta, build the vector of theta values
if prod(size(theta))==1
   theta = (0:(size(p,2)-1))* theta;
end

if length(theta) ~= size(p,2)
   error('THETA does not match the number of projections.');
end



