% Run the data propagation in frequency domain: 

% change the the data file name you want: 
datafilename = 'bub-bam1-lateral.csv';

% change to the correct file path on your computer: 
emptyfile = 'empty.csv'; 

% Propagation distance: 
PropDistance = [0.95]; % a number or vector of distance of propagation in meters

% a list of frequency at which you want to display the result: 
DisplayedFreqList = [7.0e9]; % in Hz. 



% load the incident wave (the data without objects):
[freq,incwave] = load_data_freq(emptyfile,100,50,51);  % 100: number of frequencies, 50: number of grid points in x, 51: number of grid point in y.
doDatapropagationFrequency([],datafilename,incwave,PropDistance,DisplayedFreqList)