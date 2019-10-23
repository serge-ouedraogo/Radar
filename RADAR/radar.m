clear all
clc

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m             
% Range Resolution = 1 m      
% Max Velocity = 100 m/s    
% Speed of light = 3e8
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
 
R = 110;
v = -20;

fprintf('Ground Truth Range = %d, Velocity = %d\n', R, v);

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

c = 3e8;                %speed of light in meters/sec
delta_r = 1;            %range resolution in meters
Bsweep = c/2*delta_r;   %Bsweep calculation

range_max = 200;            %given radar's max range
Tchirp = 5.5*(range_max*2/c);   %5.5 times of the trip time for maximun range

slope = Bsweep/Tchirp;          %Slope of the FMCW chirp
fprintf('Slope = %d\n', slope);

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                         
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd = 128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t = linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); %transmitted signal
Rx = zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)            
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v * t(i);
    td(i) = 2*r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+0.5*slope*t(i)^2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+0.5*slope*(t(i)-td(i))^2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT
L = Tchirp*Bsweep;  % Length of signal

% *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix=reshape(Mix,[Nr,Nd]);

% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and normalize.
signal_fft = fft(Mix,Nr);

% *%TODO* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft/L);

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:L/2+1);

% *%TODO* :
% plot FFT output 
f = Bsweep*(0:(L/2))/L;
figure('Name','FFT Output')
plot(f,signal_fft)

% Calculate range estimation
R = (c*Tchirp*f)/(2*Bsweep);

%plotting the range
figure('Name','Range from First FFT')
plot(R,signal_fft)
axis([0 200 0 0.5]);


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 12;
Td = 6;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 6;
Gd = 3;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 8;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
   % a vector to  hold signal after thresholding
   signal = zeros(Nr/2,Nd);
   cells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);

for i = 1:(Nr/2 - (2*Gr+2*Tr))
    for j = 1:(Nd - (2*Gd+2*Td))
        % Add noise in each cell  
        sum1 = sum(db2pow(RDM(i:i+2*Tr+2*Gr, j:j+2*Td+2*Gd)),'all');
        sum2 = sum(db2pow(RDM(i+Tr:i+Tr+2*Gr, j+Td:j+Td+2*Gd)),'all');    
        noise_level = sum1 - sum2;
              
        % To determine the threshold, we take the average of summed noise
        % and multiply it with the offset        
        
        threshold = noise_level/cells;      
        threshold = pow2db(threshold) + offset;
        threshold = db2pow(threshold);
        
        % Now pick the cell under test which is the T+G cells away from the
        % first training cell and seasure the signal level
        cut = db2pow(RDM(i+Tr+Gr, j+Td+Gd));
        
        % If the signal level at Cell Under Test belows the threshold then
        % assign it to 0 value
        if (cut <= threshold)
            cut = 0;
        else 
            cut = 1;
        end
        
        signal(i+Tr+Gr,j+Td+Gd) = cut;        
    end
end

   

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal);
colorbar;

