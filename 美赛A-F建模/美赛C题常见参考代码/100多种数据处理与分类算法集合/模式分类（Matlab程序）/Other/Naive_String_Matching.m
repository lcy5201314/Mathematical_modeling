function sh = Naive_String_Matching(text, x)

% Naive string matching: Find the location of the string x in text 'text'
%
% Inputs:
%	text				- Text vector
%	x					- Search string
%
% Output:
%	sh					- Locations of the string

n	= length(text);
m	= length(x);
s	= 0;
sh	= [];

while (s <= n - m),
   if strcmp(x, text(s+1:s+m)),
      disp(['The pattern occurs at shift ' num2str(s+1)])
		sh(end+1) = s+1;
   end
   s = s + 1;
end
