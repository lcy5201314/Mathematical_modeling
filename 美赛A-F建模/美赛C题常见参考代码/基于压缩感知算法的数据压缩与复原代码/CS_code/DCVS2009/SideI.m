function [ img ] = SideI( img1, img2, method)
% method = 0 -> frame interpolation
% method = 1 -> forword motion estimation
block = 16;   %block size
sr = 7;         %search range
[ height width ] = size(img2);
img = zeros( height, width);
if(method>1)
    imgb = zeros(height, width);
end
[ h1 w1 ] = size(img1);
if(h1 ~= height || w1 ~= width)
    disp('images size are not equal')
else
    switch method
        case 0  % 0 -> frame interpolation
            img = (img1+img2)/2;
        case 1  % 1 -> forword motion estimation
            [ height width ] = size(img2);
            for h=1:block:height
                for w=1:block:width
                    maxSAD = 65535;
                    for y=-sr:sr
                        for x=-sr:sr                        
                            if ( h+y>=1 && h+y<=height-block && w+x>=1 && w+x <=width-block)
                                SAD = 0;
                                for j=0:15
                                    for i=0:15
                                        SAD = SAD + abs( img2(h+j,w+i) - img1(h+y+j,w+x+i));
                                    end
                                end
                                if( SAD < maxSAD)
                                    maxSAD = SAD;
                                    mvx = x;
                                    mvy = y;
                                end
                            end
                        end
                    end
                    for j=0:15
                        for i=0:15
                            img( h+j, w+i ) = img1( h+floor(mvy/2)+j , w+floor(mvx/2)+i );
                        end
                    end
                end
            end
        case 2           %(forward MC + backward MC) /2
            [ height width ] = size(img2);
            for h=1:block:height
                for w=1:block:width
                    maxSADf  = 65535;
                    maxSADb = 65535;
                    for y=-sr:sr
                        for x=-sr:sr                        
                            if ( h+y>=1 && h+y<=height-block && w+x>=1 && w+x <=width-block)
                                SADf  = 0;
                                SADb = 0;
                                for j=0:15
                                    for i=0:15
                                        SADf  = SADf  + abs( img2(h+j,w+i) - img1(h+y+j,w+x+i));
                                        SADb = SADb + abs( img1(h+j,w+i) - img2(h+y+j,w+x+i));
                                    end
                                end
                                if( SADf < maxSADf)
                                    maxSADf = SADf;
                                    mvfx = x;
                                    mvfy = y;                                
                                end
                                if( SADb < maxSADb)
                                    maxSADb = SADb;
                                    mvbx = x;
                                    mvby = y;
                                end
                            end
                        end
                    end
                    for j=0:15
                        for i=0:15
                            img( h+j, w+i )    = img1( h+floor(mvfy/2)+j , w+floor(mvfx/2)+i );
                            imgb( h+j, w+i ) = img2( h+floor(mvby/2)+j , w+floor(mvbx/2)+i );
                        end
                    end
                end
            end
            img = (img + imgb)/2;
    end
end
