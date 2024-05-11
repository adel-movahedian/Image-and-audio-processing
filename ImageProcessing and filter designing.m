%% Q1 match filter
%% 1.2 
clear; clc; close all;
Slength = input("what is Slength?");
Nlength = input("what is Nlength?");
signal = zeros(1, Slength); %genereating signal
symbol = complex(randn(1,Nlength) , randn(1,Nlength));%unsing complex func for generating symbol(converting real and imaginary part to a complex number)
noise = 0.5 * complex(randn(1, Slength) , randn(1, Slength));%noise amplitude should be half of the symbol amplitude
segment_start = randi([1, Slength - Nlength + 1], 1, 2); %generating signal with symbol in random specific places
for i = 1:Nlength
    signal(segment_start(1) + i - 1) = symbol(i);
    signal(segment_start(2) + i - 1) = symbol(i);
end
sample = signal + noise ;   %sample would be the sum of signal and nise

%showing signals :
f1 = figure('Name','Generated Symbol','NumberTitle','off','Position', [100 350 560 220]);
plot(1:Nlength , abs(symbol) , 'black');
title('Generated Symbol');
grid minor
f2 = figure('Name','Generated Sample','NumberTitle','off','Position', [100 50 560 220]);
plot(1:length(sample) , abs(sample) , 'black');
title('Generated Sample');
grid minor
% 1.3
flipedSymbol = flip(symbol) ; %reversing signal
matchedFilter = conj(flipedSymbol) ; %The conj function is used to take the complex conjugate of the symbol array
filteredSignal = conv(signal , matchedFilter , 'same') ;%The 'same' argument tells conv to return a central portion of the convolution that is the same size as the signal array. This ensures that the output has the same length as the input
normalizedfilteredSignal = filteredSignal / norm(filteredSignal) ;%normalizing the signal
f3 = figure('Name','filtered signal','NumberTitle','off','Position', [700 100 560 440]);
plot(1:length(normalizedfilteredSignal) , abs(normalizedfilteredSignal) , 'red' ) ;
title('Filtered signal');
grid minor
indices = find(normalizedfilteredSignal > 0.3);
disp("symbol is on : "+ indices);
%% 1.4

clear; clc; close all;
%Loading samples and symbols
load('sample1.mat');
load('sample2.mat');
load('sample3.mat');
load('sample4.mat');
load('sample5.mat');
load('sample6.mat');
load('sample7.mat');
load('sample8.mat');
load('symbol1.mat');
load('symbol2.mat');

% sample1 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample1 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filteredsignal(sample1&symbol1)');

% sample1 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample1 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample1 & symbol2)');

% sample2 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample2 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample2 & symbol1)');

% sample2 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample2 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample2 & symbol2)');

% sample3 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample3 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample3 & symbol1)');

% sample3 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample3 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample3 & symbol2)');

% sample4 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample4 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample4 & symbol1)');

% sample4 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample4 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample4 & symbol2)');

% sample5 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample5 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample5 & symbol1)');

% sample5 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample5 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample5 & symbol2)');

% sample6 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample6 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (ample6 & symbol1)');

% sample6 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample6 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample6 & symbol2)');

% sample7 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample7 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample7 & symbol1)');

% sample7 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample7 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample7 & symbol2)');

% sample8 & symbol1
flipedSymbol = flip(symbol1) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample8 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample8 & symbol1)');

% sample8 & symbol2
flipedSymbol = flip(symbol2) ;
matchedFilter = conj(flipedSymbol) ;
filteredSignal = conv(sample8 , matchedFilter , 'same') ;
filteredSignal = filteredSignal / norm(filteredSignal) ;
figure
plot(1:length(filteredSignal) , abs(filteredSignal) , 'b' ) ;
title('Filtered signal (sample8 & symbol2)');

%% Q2 convolution
 %% 2.2.1 
%the func is at the end of the code

 %% 2.2.2
clc; clear; close all;
%signals forms :
f1 = figure('Name','Signals','NumberTitle','off');
f1.Position (1:2) = [0 200] ;

subplot (3 , 1 , 1);
n = -20 : 1 : 20;
u1 = n >= 6;
u2 = n >= 0;
x1 =u2 - 4*n.*u1 + (n==-5) - (n==-6);
stem (n , x1 ,'Color','red');
title ('$x_1[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 2);
n = -20 : 1 : 20;
u1 = n >= 0;
u2 = n >= 5;
x2 = (0.5) .^ n .* (u1 - u2);
stem (n , x2 , 'Color','blue');
title ('$x_2[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 3);
n = -20 : 1 : 20;
u1 = n >= 0;
x3 = exp(sin((pi/5) .* n))  .* (u1 );
stem (n , x3 , 'Color','black');
title ('$x_3[n]$','Interpreter','latex');
grid minor;

% convolve : ( myConv )

f2 = figure('Name','Convolutions','NumberTitle','off');
f2.Position (1:2) = [600 200] ;

subplot (3 , 1 , 1);
n = -40 : 1 : 40;
x4 = myConv (x1 , x2);
stem (n , x4 , 'Color','green');
title ('$x_1[n] * x_2[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 2);
n = -40 : 1 : 40;
x5 = myConv (x1 , x3);
stem (n , x5 , 'Color',"#77AC30");
title ('$x_1[n] * x_3[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 3);
n = -40 : 1 : 40;
x6 = myConv (x2 , x3);
stem (n , x6 , 'Color','m');
title ('$x_2[n] * x_3[n]$','Interpreter','latex');
grid minor;
%% 2.2.3 
clc; clear; close all;
%signals forms :
n = -20 : 1 : 20;
u1 = n >= 6;
u2 = n >= 0;
x1 =u2 - 4*n.*u1 + (n==-5) - (n==-6);

n = -20 : 1 : 20;
u1 = n >= 0;
u2 = n >= 5;
x2 = (0.5) .^ n .* (u1 - u2);

n = -20 : 1 : 20;
u1 = n >= 0;
x3 = exp(sin((pi/5) .* n))  .* (u1 );

% convolve : ( conv )

f2 = figure('Name','Convolutions','NumberTitle','off');
f2.Position (1:2) = [600 200] ;

subplot (3 , 1 , 1);
n = -40 : 1 : 40;
x4 = conv (x1 , x2);
stem (n , x4 , 'Color','green');
title ('$x_1[n] * x_2[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 2);
n = -40 : 1 : 40;
x5 = conv (x1 , x3);
stem (n , x5 , 'Color',"#77AC30");
title ('$x_1[n] * x_3[n]$','Interpreter','latex');
grid minor;

subplot (3 , 1 , 3);
n = -40 : 1 : 40;
x6 = conv (x2 , x3);
stem (n , x6 , 'Color','m');
title ('$x_2[n] * x_3[n]$','Interpreter','latex');
grid minor;
figure;
 %second part 
n = -20 : 1 : 20;
a = randi( [0 3] , 1 , 21); %generating random arrays with desired sizes
b = randi( [0 3] , 1 , 21);

tic;
subplot (2 , 1 , 1);
stem (n , conv( a , b) , 'Color','m');
title ('$Using\  conv$','Interpreter','latex');
grid minor;
toc;

tic;
subplot (2 , 1 , 2);
stem (n , myConv( a , b) , 'Color','r');
title ('$Using\  myConv$','Interpreter','latex');
grid minor;
toc;
%% 2.3
clc; clear; close all;
A = [1 1 0 ; 1 0 1 ; 0 1 1]
B = [1 2 3 4 5 ; 6 7 8 9 10 ; 11 12 13 14 15 ; 16 17 18 19 20 ; 21 22 23 24 25]
Convolution = conv2 (A , B)
%% Image processing
%% 3.1 Image Processing - Some Basic Commands
%1
clc; clear; close all;
f1 = figure('Name','RGB','NumberTitle','off');
f1.Position (1 , 1);

X = imread ('flower.jpg');
black = zeros ( size ( X , [1,2] ) );

subplot (2 , 2 , 1);
image (X);

subplot (2 , 2 , 2);
R = cat (3 , X( : , : , 1) , black , black );
image(R); 

subplot (2 , 2 , 3);
G = cat (3 , black , X( : , : , 1) , black );
image(G); 

subplot (2 , 2 , 4);
B = cat (3 , black , black , X( : , : , 3) );
image(B); 

%2

f2 = figure('Name','figure2','NumberTitle','off');
mean = ( X( : , : , 1) + X( : , : , 2) + X( : , : , 3) ) / 3;
imshow (mean);
title ('$mean$','Interpreter','latex');
f3 =  figure('Name','figure3','NumberTitle','off');
gray = rgb2gray (X);
imshow (gray);
title ('$rg2gray$','Interpreter','latex');

%3
figure;
d = im2double (mean);
imshow (d);
title ('$im2double$','Interpreter','latex');
%% Image Processing - Edge Detection (B&W)

clc; clear; close all;

f1 = figure('Name','Original Picture','NumberTitle','off');
f1.Position (1:2) = [-250 500] ;

X = imread ('flower.jpg');
IMG = im2double ( rgb2gray (X) );
imshow (X);
title ('$Original\  Picture$','Interpreter','latex');

f2 = figure('Name','Horizontal Edges','NumberTitle','off');
f2.Position (1:2) = [250 500] ;

EDKx = [ -1 0 1 ; -2 0 2 ; -1 0 1 ];
IMGx = conv2 (IMG , EDKx);
imshow (IMGx);
title ('$Horizontal\  Edges$','Interpreter','latex');

f3 = figure('Name','Vertical Edges','NumberTitle','off');
f3.Position (1:2) = [530 500] ;

EDKy = [ -1 -2 -1 ; 0 0 0 ; 1 2 1 ];
IMGy = conv2 (IMG , EDKy);
imshow (IMGy);
title ('$Vertical\  Edges$','Interpreter','latex')

f4 = figure('Name','All Edges','NumberTitle','off');
f4.Position (1:2) = [950 500] ;

IMGxy = sqrt (IMGx .^ 2 + IMGy .^ 2);
imshow (IMGxy);
title ('$All\  Edges$','Interpreter','latex')
%% Image Processing - Edge Detection 

clc; clear; close all;

X = imread ('flower.jpg');
EDKx = [ -1 0 1 ; -2 0 2 ; -1 0 1 ];
EDKy = [ -1 -2 -1 ; 0 0 0 ; 1 2 1];

%dividing channels
r = im2double ( X ( : , : , 1) );
g = im2double ( X ( : , : , 2) );
b = im2double ( X ( : , : , 3) );

%horizontal
IMGrh = conv2 ( EDKx , r);
IMGgh = conv2 ( EDKx , g);
IMGbh = conv2 ( EDKx , b);

%vertical
IMGrv = conv2 ( EDKy , r);
IMGgv = conv2 ( EDKy , g);
IMGbv = conv2 ( EDKy , b);

IMGr = sqrt (IMGrh .^ 2 + IMGrv .^ 2);
IMGg = sqrt (IMGgh .^ 2 + IMGgv .^ 2);
IMGb = sqrt (IMGbh .^ 2 + IMGbv .^ 2);

IMG = cat (3 , IMGr , IMGg , IMGb);

imshow (IMG);
%% functions
function output = myConv (a , b)
    % These 6 lines are for preparing signals for convolution operation 
    a = flip (a);
    x = length(b) : length(b)+length(a)-1;
    c = zeros (1 , length(a) + 2 * (length(b) - 1) );
    c(x) = a(x - (length(b)-1));
    d (1 : length(b)) = b ;
    d ( length(a) + 2 * (length(b) - 1) ) = 0;
    
    for i = 1 : length(a)+length(b)-1 
        g = c .* d;
        output(i) = sum(g);

        if i<(length(a)+length(b)-1)
            % shifting our signal to right
            s = i+1 : length(d);
            d(s) = d(s-1);
            d(1:i) = 0;
        end
    end
    output = flip (output);
end

