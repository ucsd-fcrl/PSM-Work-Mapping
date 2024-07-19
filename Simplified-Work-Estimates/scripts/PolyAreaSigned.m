function area = PolyAreaSigned(x,y,dim)
%POLYAREA Area of polygon.
%   POLYAREA(X,Y) returns the area of the polygon specified by
%   the vertices in the vectors X and Y.  If X and Y are matrices
%   of the same size, then POLYAREA returns the area of
%   polygons defined by the columns X and Y.  If X and Y are
%   arrays, POLYAREA returns the area of the polygons in the
%   first non-singleton dimension of X and Y.
%
%   The polygon edges must not intersect.  If they do, POLYAREA
%   returns the difference between the clockwise
%   encircled areas and the counterclockwise encircled areas.
%
%   POLYAREA(X,Y,DIM) returns the area of the polygons specified
%   by the vertices in the dimension DIM.
%
%   Class support for inputs X,Y:
%      float: double, single

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.3 $  $Date: 2010/08/23 23:11:57 $

if nargin==1
	error(message('MATLAB:polyarea:NotEnoughInputs'));
end

if ~isequal(size(x),size(y))
	error(message('MATLAB:polyarea:XYSizeMismatch'));
end

if nargin==2,
	[y,nshifts] = shiftdim(y);
	x = shiftdim(x);
elseif nargin==3,
	perm = [dim:max(length(size(y)),dim) 1:dim-1];
	x = permute(x,perm);
	y = permute(y,perm);
end

siz = size(y);
if ~isempty(y),
	area = reshape(sum( (y([2:siz(1) 1],:) - y(:,:)).* ...
		(x([2:siz(1) 1],:) + x(:,:)))/2,[1 siz(2:end)]);
else
	area = sum(y); % SUM produces the right value for all empty cases
end

if nargin==2,
	area = shiftdim(area,-nshifts);
elseif nargin==3,
	area = ipermute(area,perm);
end

% PolyAreaSigned
%
% Computes the signed area of a polygon specified by a set of ordered
% points. If the polygon in self intersecting returns the difference of
% area of positive and negative loops.

% function polyArea= PolyAreaSigned(x,y);
% 	xPermute = [x(2:end) x(1)];
% 	polyArea = -0.5*sum((x + xPermute) .* [diff(y) (y(1)-y(end))]);
