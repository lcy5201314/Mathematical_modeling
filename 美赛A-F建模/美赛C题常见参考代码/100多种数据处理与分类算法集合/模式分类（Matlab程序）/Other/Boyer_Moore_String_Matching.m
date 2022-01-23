function sh = Boyer_Moore_String_Matching(text, x)

% Boyer-Moore string matching: Find the location of the string x in text 'text'
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
   j = m;
   while ((j > 0) & (x(j) == text(s+j))),
      j = j - 1;
   end
   if (j == 0),
      disp(['The pattern occurs at shift ' num2str(s+1)])
      sh(end+1) 	= s+1;
      s				= s + G(x, 0);
   else
      s				= s + max(G(x, j), j - F(x, text(s + j)));
   end
end
%END

function s = G(x, m)
%Good suffix function: Find the leftmost occurence of the 'm' characters
n = length(x);

s = 0;
for i = n-2*m+1:-1:1,
   if strcmp(x(end-m+1:end), x(i:i+m-1)),
      s = i;
      break
   end
end

function s = F(x, c)
%Last occurrence function
n = length(x);

s = 0;
for i = n:-1:1,
   if (x(i) == c),
      s = i;
   end
end
