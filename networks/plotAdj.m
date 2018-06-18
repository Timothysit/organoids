function plotAdj(adjM, goodElectrodes)
% note that this is quite customised for the mecp2 project (MEA 8 x 8) data
% INPUT 
    % plotAdj | N x N adjacency matrix (where N is the number of active electrodes) 
    % goodElectrodes | N x 1 vector saying which electrodes are active;
    % only plot those

% Network adjacency plot, like the one by Manuel 
% author: Tim Sit github/timothysit 
% last update: 20180407 
% TODO: check that the coordinates agrees with other heatmap visualisations

%% generate the possible 2-val combinations of 1 to 8, with repetition

% n = 8; k = 2;
% nk = nmultichoosek(1:n,k);
% coord=zeros(0,k);
% for i=1:size(nk,1)
%     pi = perms(nk(i,:));
%     coord = unique([coord; pi],'rows');
% end
% 
% coord(:, 1) = 9 - coord(:, 1); 

% coord needs to count from top to bottom columnwise, ie. [8, 1], [8,
% 2]... 

yTemp = 1:8; 
yCoord = repmat(fliplr(yTemp), 1, 8);

xTemp = 1:8; 
xCoord = repelem(xTemp, 8);

coord = [xCoord', yCoord'];

%% prune coord 


coord(coord(:, 1) == 1 & coord(:, 2) == 1, :) = [];
coord(coord(:, 1) == 1 & coord(:, 2) == 8, :) = [];
coord(coord(:, 1) == 8 & coord(:, 2) == 8, :) = [];
coord(coord(:, 1) == 8 & coord(:, 2) == 1, :) = [];

coord = coord(goodElectrodes, :);


%% do some pruning of the adjacency so we don't plot everytying 

edges = adjM; % flip things to make it agree with the controllabiliyt heatmap
% edges(edges == 1) = 0; % remove self-correlation 
% threshold = prctile(edges(:), 75); % get the 95th percentile 
threshold = 0.75;
edges(edges < threshold) = 0.0001; 
edges(isnan(edges)) = 0.0001; % I think 0 doens't quite work

%% Custom colormap 

% generated from here: 
% http://jdherman.github.io/colormap/

C = [255,255,255;
253,255,252;
252,254,250;
250,254,247;
248,253,245;
247,253,242;
245,252,240;
243,252,237;
242,251,234;
240,251,232;
238,250,229;
237,250,227;
235,249,224;
233,249,222;
232,248,219;
230,248,216;
228,247,214;
227,247,211;
225,246,209;
223,246,206;
222,245,204;
220,245,201;
218,244,198;
217,244,196;
215,243,193;
213,243,191;
212,242,188;
210,242,186;
208,241,183;
206,241,180;
205,240,178;
203,240,175;
201,239,173;
200,239,170;
198,239,168;
196,238,165;
195,238,162;
193,237,160;
191,237,157;
190,236,155;
188,236,152;
186,235,150;
185,235,147;
183,234,145;
181,234,142;
180,233,139;
178,233,137;
176,232,134;
175,232,132;
173,231,129;
171,231,127;
170,230,124;
168,230,121;
166,229,119;
165,229,116;
163,228,114;
161,228,111;
160,227,109;
158,227,106;
156,226,103;
155,226,101;
153,225,98;
151,225,96;
150,224,93;
148,223,92;
147,220,94;
145,216,96;
144,213,98;
143,209,100;
141,206,102;
140,202,104;
139,199,106;
138,195,108;
136,192,111;
135,188,113;
134,185,115;
132,181,117;
131,178,119;
130,174,121;
128,170,123;
127,167,125;
126,163,127;
124,160,129;
123,156,132;
122,153,134;
121,149,136;
119,146,138;
118,142,140;
117,139,142;
115,135,144;
114,132,146;
113,128,148;
111,125,150;
110,121,153;
109,118,155;
107,114,157;
106,111,159;
105,107,161;
104,104,163;
102,100,165;
101,97,167;
100,93,169;
98,90,171;
97,86,173;
96,83,176;
94,79,178;
93,76,180;
92,72,182;
90,69,184;
89,65,186;
88,62,188;
87,58,190;
85,54,192;
84,51,194;
83,47,197;
81,44,199;
80,40,201;
79,37,203;
77,33,205;
76,30,207;
75,26,209;
73,23,211;
72,19,213;
71,16,215;
69,12,218;
68,9,220;
67,5,222;
66,2,224;
66,0,223;
67,0,220;
68,0,216;
70,0,213;
71,0,209;
73,0,205;
74,0,202;
76,0,198;
77,0,195;
78,0,191;
80,0,188;
81,0,184;
83,0,181;
84,0,177;
86,0,174;
87,0,170;
89,0,167;
90,0,163;
91,0,160;
93,0,156;
94,0,153;
96,0,149;
97,0,146;
99,0,142;
100,0,138;
101,0,135;
103,0,131;
104,0,128;
106,0,124;
107,0,121;
109,0,117;
110,0,114;
111,0,110;
113,0,107;
114,0,103;
116,0,100;
117,0,96;
119,0,93;
120,0,89;
121,0,86;
123,0,82;
124,0,78;
126,0,75;
127,0,71;
129,0,68;
130,0,64;
131,0,61;
133,0,57;
134,0,54;
136,0,50;
137,0,47;
139,0,43;
140,0,40;
141,0,36;
143,0,33;
144,0,29;
146,0,26;
147,0,22;
149,0,19;
150,0,15;
151,0,11;
153,0,8;
154,0,4;
156,0,1;
157,0,0;
158,0,0;
159,0,0;
160,0,0;
161,0,0;
163,0,0;
164,0,0;
165,0,0;
166,0,0;
167,0,0;
168,0,0;
169,0,0;
170,0,0;
171,0,0;
173,0,0;
174,0,0;
175,0,0;
176,0,0;
177,0,0;
178,0,0;
179,0,0;
180,0,0;
181,0,0;
183,0,0;
184,0,0;
185,0,0;
186,0,0;
187,0,0;
188,0,0;
189,0,0;
190,0,0;
191,0,0;
193,0,0;
194,0,0;
195,0,0;
196,0,0;
197,0,0;
198,0,0;
199,0,0;
200,0,0;
202,0,0;
203,0,0;
204,0,0;
205,0,0;
206,0,0;
207,0,0;
208,0,0;
209,0,0;
210,0,0;
212,0,0;
213,0,0;
214,0,0;
215,0,0;
216,0,0;
217,0,0;
218,0,0;
219,0,0;
220,0,0;
222,0,0;
223,0,0;
224,0,0;
225,0,0;
226,0,0;
227,0,0];


%% actual network plot 
% figure
% gplot(a,coord) % most simple way of doing it
% h = adjacency_plot_und(a, coord); % requires brain connectivity toolbox

weights = sum(edges) / sum(sum(edges));
wgPlot(edges, coord, 'edgeColorMap', colormap(C/255), 'vertexWeight', weights', 'vertexScale',200,'edgeWidth',2);
% caxis([0 1])
% vertex scale oringal value: 200
% edgewidth original value: 2

% TS 20180612: try out different colormaps
% note the originla black and white uses: colormap(flipud(gray(256)))

aesthetics

end 