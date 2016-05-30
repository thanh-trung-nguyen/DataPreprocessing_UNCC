function u2 = extract_target_signal(u)


threshold = 0.3; % threshold for curve compared to the maximum curve
Threshold2 = 0.8; % threshold for minimum peak

[Nt,Np] = size(u);
u2 = 0*u;
Min = min(min(u));

for n = 1:Np
    v = u(:,n);
    Minn = min(min(v));
    if Minn < threshold*Min
        Peak = find(v<Threshold2*Minn,1,'first');
        while v(Peak) <0 && Peak > 1
            Peak = Peak - 1;
        end
        while v(Peak) >= 0 && Peak > 1
            Peak = Peak - 1;
        end
        while v(Peak) < 0 && Peak > 1
            Peak = Peak - 1;
        end
        v(1:Peak) = 0;
        u2(:,n) = v;
    end    
end
% cut out the positive peaks at the beginning:
for n = 1:Np
    k = 1;
    while (k < Nt-2) && u2(k,n) >=0       
        u2(k,n) = 0;
        k = k+1;
    end
end