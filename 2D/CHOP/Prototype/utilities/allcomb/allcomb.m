function A = allcomb(varargin)

% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
%    of the number of elements of the N inputs. This functionality is also
%    known as the Cartesian Product. The arguments can be numerical and/or
%    characters, or they can be cell arrays.
%
%    Examples:
%       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
%       % -> [ 1  -3   0
%       %      1  -3   1
%       %      1   8   0
%       %        ...
%       %      5  -3   1
%       %      5   8   1 ] ; % a 12-by-3 array
%
%       allcomb('abc','XY') % character arrays
%       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
%
%       allcomb('xy',[65 66]) % a combination
%       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
%
%       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
%       % -> {  'hello'  'Joe'        [99999]
%       %       'hello'  'Joe'             []
%       %       'hello'  [1x3 double] [99999]
%       %       'hello'  [1x3 double]      []
%       %       'Bye'    'Joe'        [99999]
%       %       'Bye'    'Joe'             []
%       %       'Bye'    [1x3 double] [99999]
%       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
%
%    ALLCOMB(..., 'matlab') causes the first column to change fastest which
%    is consistent with matlab indexing. Example: 
%      allcomb(1:2,3:4,5:6,'matlab') 
%      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
%
%    If one of the arguments is empty, ALLCOMB returns a 0-by-N empty array.
%    
%    See also NCHOOSEK, PERMS, NDGRID
%         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)

% for Matlab R2011b
% version 4.0 (feb 2014)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
%     added option to let the first input run fastest (suggestion by JD)
% 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
%     output by NDGRID. Thanks to Jan Simon, for showing this little trick
% 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
% return an empty output.
% 2.1 (feb 2011). A cell as input argument caused the check on the last
%      argument (specifying the order) to crash.
% 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
% 3.0 (may 2012) removed check for doubles so character arrays are accepted
% 4.0 (feb 2014) added support for cell arrays

NC = nargin ;
% default: enter arguments backwards, so last one (AN) is changing fastest
ii = NC:-1:1 ;

% check for empty inputs
if NC > 1
    % non-cell input, assuming all numerical values or strings
    % flip using ii if last column is changing fastest
    [A{ii}] = ndgrid(varargin{ii}) ;
    % concatenate
    A = reshape(cat(NC+1,A{:}),[],NC) ;
elseif NC==1
    A = varargin{1}(:) ; % nothing to combine
end

