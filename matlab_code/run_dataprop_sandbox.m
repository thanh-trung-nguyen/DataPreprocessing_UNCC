% frequency data propagation for the sandbox


% load data of object and sand alone: 
datafile = '/home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-aluminum-cylinder.csv';

[freq,dat3d,dat2d] = load_data_freq(datafile,300,50,51);

datafilesand = '/home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-sand-only.csv';

[freq,sand3d,sand2d] = load_data_freq(datafilesand,300,50,51);



% propagate data:
LightSpeed = 299792458; % m/s.

PropDistance = 0.57; % propagate to the front surface of the sand box.
x = linspace(-0.5,0.5,50); 
y = linspace(-0.5, 0.5, 51); 

datprop = wave_propagation_freq_3d(dat3d-sand3d,PropDistance,freq',x,y,LightSpeed); 


for n = 1:200; 
    imagesc(abs(reshape(datprop(n,:),50,51))'); 
    title(freq(n)); 
    pause(0.2); 
end