% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons

%% SET UP THE NETWORK

% 1- NETWORK SIZE:
Ne=800; Ni=200; % Excitatory, inhibitory. Ne+Ni is total neurons.

% 2 - GLOBAL PARAMETERS THAT SET OUR NEURON MODEL. DEFAULT IS SPIKING
% NEURON:
% Set initial conditions of neurons, with some variability provided by the
% vectors re and ri containing random numbers between 0 and 1.
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];

% 3 - SET UP THE CONNECTIVITY MATRIX: DIRECTED NETWORK
% In this construction, 1=connection exist, 0=no connection.
% Connectivity is set as random. Then, a fraction of connections are set 0.
% Note that effectively this is an Erdös-Rényi graphs, with no spatial
% characteristics. % If you want to symmetrize the network to make 
% undirected networks (not realistic in neuroscience), use A = (A + A.')/2
% after line 27.
frac_delete=0.6; % Set this fraction of connections to zero.
A=[rand(Ne+Ni)];
A(A<frac_delete)=0;
A(A>0)=1;
A = A - diag(diag(A)); % Make the diagonal elements 0 (no self-connection).
%figure;
%imagesc(A);

% 4 - SET SYNAPTIC WEIGHTS (STRENGTHS) OF CONNECTIONS.
% I represents the amnplitudes of evoked currents (EPSC and IPSC). 
% Increase the max. excitatory weigth to induce synchronization.
% Be careful. The balance between excitation and inhibition governs the
% dynamics of the network!
MAX_EXC_WEIGTH=1; MAX_INH_WEIGTH=.5; 
W=[MAX_EXC_WEIGTH*rand(Ne+Ni,Ne), -MAX_INH_WEIGTH*rand(Ne+Ni,Ni)];

% 5 - The final connectivity matrix S is the element-wise multiplication of A
% and W.
S=A.*W;  % Note that S is directed and weighted!!
S0 = A.*W; %initial S
% 6 - DEFINE NOISE STRENGTH. Remember that noise (or random inputs) is the
% main drive of spontaneous activity. Default is around 5. 
NOISE_MAX=2; %Initial value in code = 5 


%% MAIN SIMULATION:

v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u
firings=[]; % spike timings
percent = []

%CODI NOSTRE
last_spike_times = zeros(Ne+Ni,1);
beta = 1 %cnt damping coefficient
thau = 200 %cnt recovery time
total = [];
for t=1:1000 % simulation of 1000 ms
I=[NOISE_MAX*randn(Ne,1);2*randn(Ni,1)]; % NOISE or thalamic input
fired=find(v>=30); % indices of spikes
%S(fired,:) = 
firings=[firings; t+0*fired,fired];
v(fired)=c(fired);
u(fired)=u(fired)+d(fired);

%CODI NOSTRE
%Per a aquest temps que ha fet fire mirem, de la neurona, l'últim firing al
%vector firings. Si el interval de temps és inferior a thau, multipliquem
%per beta
%t - t_ub < 10 multipliquem per beta
restored = find((t+0*last_spike_times)-last_spike_times >= thau);
reduced = find((t+0*fired)-last_spike_times(fired) <= thau);
S(restored, :) = S0(restored,:);
S(reduced, :) = beta*S(reduced, :);
test = sum(sum(S));
total = [total;test]
I=I+sum(S(:,fired),2);%sums all the noise plus the inputs from the neighbouring neurons (vector of size of v)
v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
u=u+a.*(b.*v-u); % stability

last_spike_times(fired) = t; %serveix per a pròxim pas

%Grafiquem quin % fa firing
percent = [percent, length(fired)/(Ne+Ni)];


end

figure;
time_vector = [1:1000]
plot(transpose(time_vector), total)
%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
figure;
scatter(firings(:,1),firings(:,2)); 
xlabel('time (ms)');
ylabel('neuron');

time = linspace(1,1000, 1000)
% Percent plot. Time on X, percent on Y.
figure;
plot(time, percent)
xlabel('time (ms)');
ylabel('percent');