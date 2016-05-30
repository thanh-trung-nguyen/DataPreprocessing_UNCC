function FirstSampleIdx = find_first_sample(v,Threshold2)
% find the first arriving sample of the target's signal at a given
% time-dependent curve

Minn = min(v);
Peak = find(v<Threshold2*Minn,1,'first');
while v(Peak) < 0 && Peak > 1
    Peak = Peak - 1;
end

while v(Peak) >= 0 && Peak > 1
    Peak = Peak - 1;
end
while v(Peak) < 0 && Peak > 1
    Peak = Peak - 1;
end
FirstSampleIdx = Peak;
