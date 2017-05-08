function order  = survey_transectorder(pos)

        

% optimize the transect order

% Get the positions
n = length(pos.pos);
lats = inf(2*n+1,1);
lons = inf(2*n+1,1);
r = inf(2*n+1);

%% Organize the positions
lats(1) = pos.startpos.lat;
lons(1) = pos.startpos.lon;
for i=1:n
    lats([2*i 2*i+1]) = [pos.pos(i).latstart pos.pos(i).latstop];
    lons([2*i 2*i+1]) = [pos.pos(i).lonstart pos.pos(i).lonstop];
    r_t(i) = pos.pos(i).R(end);%nmi
end

%% calculate the distances between transects
for i=1:size(r,1)
    for j=i:size(r,2)
            r(i,j) = m_lldist([lons(i) lons(j)],[lats(i) lats(j)])./1852;%km
    end
end
%% combine transcts

%transect order
order = NaN(n,1);
order(1)=1;
notused = true(size(r,1),1);
notused(1)=false;

for i=1:n
%    disp(notused')
    % find closest point
    [~,ind] = min(r(i,notused));
    % Ind related to the original vector
    dum1=find(cumsum(notused)==ind);
    ind2=dum1(1);
    
    % And find other side of transect
    if mod(ind2,2)
        ind2(2)=ind2(1)-1;
        order(i)=-(min(ind2)/2);
    else
        ind2(2)=ind2(1)+1;
        order(i)=(min(ind2)/2);
    end
    notused(ind2)=false;
end


% 
switch pos.vessel
    case 'NO'
        order=[1 -2 -12 -14 15 13 3 4 5 6 -7 -8 -16 -17 -9 -10 -18 -19 -11];
    case 'RU'
        order = 1:16;
    case 'IS'
        order = [-2 -3 -4 -5 6 7 1 8  -9 -10 -11];
    case 'EU'
        order = [-1 -2 -3  -4  -5  -6  -8 7];
    case 'FO'
        order = [1 2 3 4 -5 -6];
end
