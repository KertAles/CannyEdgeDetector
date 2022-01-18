%%
% images | derivatives | ang_mag | nonmax | LT_HT | hysteresis
picturesque = [false false false false false false];

I_raw = imread('./images/0001.png');

I = imgaussfilt(medfilt2(I_raw), 0.5);

if picturesque(1)
    figure()
    tiledlayout(1,2)

    nexttile
    imshow(I_raw)
    title('Raw image')

    nexttile
    imshow(I)
    title('Smooth')
end
%%
filter = [-1 0 1
          -2 0 2
          -1 0 1];

%filter = [-1 1];

I_x = conv2(I, filter, 'same');
I_y = conv2(I, filter', 'same');

if picturesque(2)
    figure()
    tiledlayout(1,2)
   
    nexttile
    imshow(uint8(I_x))
    title('X der')
    
    nexttile
    imshow(uint8(I_y))
    title('Y der')
end

%%

I_grad = atan2(I_y, I_x);
I_mag = sqrt(I_x.^2 + I_y.^2);

if picturesque(3)
    figure()
    tiledlayout(1,2)

    nexttile
    imshow(I_grad)
    title('Gradient Angle')

    nexttile
    imshow(uint8(I_mag))
    title('Gradient Magnitude (Pop pop)')
end

%%

I_filtered = zeros(size(I));

for i = 2:size(I, 1)-1
    for j = 2:size(I, 2)-1
        ang_type = get_edge_type(I_grad(i, j));
        val1 = 0;
        val2 = 0;

        if ang_type == 0 
            val1 = I_mag(i, j+1);
            val2 = I_mag(i, j-1);
        elseif ang_type == 1
            val1 = I_mag(i-1, j-1);
            val2 = I_mag(i+1, j+1);
        elseif ang_type == 2
            val1 = I_mag(i+1, j);
            val2 = I_mag(i-1, j);
        elseif ang_type == 3
            val1 = I_mag(i-1, j+1);
            val2 = I_mag(i+1, j-1);
        end

        if val1 < I_mag(i, j) && val2 < I_mag(i, j)
            I_filtered(i, j) = I_mag(i, j);
        end
    end
end


if picturesque(4)
    figure()
    tiledlayout(1,2)

    nexttile
    imshow(uint8(I_mag))
    title('Raw magnitude')

    nexttile
    imshow(uint8(I_filtered))
    title('Filtered magnitude(nonmaxima)')
end


%%

%otsu = graythresh(I_filtered);
%max_val = max(max(I_filtered));

%LT = otsu / 10 * max_val;
%HT = otsu / 5 * max_val;
%%
LT = median(median(I_filtered(I_filtered > 0))) * 1.1;
HT = LT * 3;

mag_LT = (I_filtered>=LT);
mag_HT = (I_filtered>=HT);

if picturesque(5)
    figure()
    tiledlayout(1,2)

    nexttile
    imshow(mag_LT)
    title('Low thresh')

    nexttile
    imshow(mag_HT)
    title('High thresh')
end



mag_hyst = mag_HT;
nu_mag_LT = mag_LT - mag_HT;


for i = 2:size(I, 1)-1
    for j = 2:size(I, 2)-1
        if mag_HT(i,j) > 0

            for k = i-1:i+1
                for l = j-1:j+1
                    if nu_mag_LT(k, l) > 0 
                        mag_hyst(k, l) = 1;
                        [nu_mag_LT, mag_hyst] = do_hysteresis(mag_hyst, nu_mag_LT, mag_HT, k, l);
                    end
                end
            end

        end
    end
end


if picturesque(6)
    figure()
    tiledlayout(1,3)

    nexttile
    imshow(uint8(I_raw))
    title('Original')

    nexttile
    imshow(uint8(I_filtered))
    title('Filtered')

    nexttile
    imshow(uint8(mag_hyst))
    title('Hystersised')
end
%%
mat_edge = edge(I,"canny_old");

if true
    figure()
    tiledlayout(1,3)

    nexttile
    imshow(uint8(I_raw))
    title('Original')

    nexttile
    imshow(mag_hyst)
    title('Our result')

    nexttile
    imshow(mat_edge)
    title('Matlab Edge')
end

mat_logic = mat_edge > 0;

TP = sum(mat_logic & mag_hyst);
FP = sum(mag_hyst) - TP;
FN = sum(~mag_hyst & mat_logic);
TN = sum(~mag_hyst) - FN;

acc = (TP + TN) / (TP + FP + FN + TN);
pr = TP / (TP + FP);
re = TP / (TP + FN);
f1 = 2 / (1/pr + 1/re);


%%

function [type] = get_edge_type(angle)
    type = -1;
    % Used types :
    % 0 - horizontal, 1 - +45, 2 - vertical, 3 - -45
    if (angle <= pi/8 && angle >= -1*pi/8) || (angle >= 7*pi/8 || angle <= -7*pi/8)
        type = 0;
    elseif (angle >= pi/8 && angle <= 3*pi/8) || (angle >= -7*pi/8 && angle <= -5*pi/8)
        type = 1;
    elseif angle >= 3*pi/8 && angle <= 5*pi/8 || (angle >= -5*pi/8 && angle <= -3*pi/8)
        type = 2;
    elseif angle >= 5*pi/8 && angle <= 7*pi/8 || (angle >= -3*pi/8 && angle <= -1*pi/8)
        type = 3;
    end

end



function [mag_LT, mag_hyst] = do_hysteresis(mag_hyst, mag_LT, mag_HT, i, j)
    mag_LT(i, j) = 0;
    
    for k = i-1:i+1
        for l = j-1:j+1
            if k > 0 && l > 0 && k <= size(mag_HT, 1) && l <= size(mag_HT, 2)
                if mag_LT(k, l) > 0
                    mag_hyst(k, l) = 1;
                    [mag_LT, mag_hyst] = do_hysteresis(mag_hyst, mag_LT, mag_HT, k, l);
                end
            end
        end
    end
end


