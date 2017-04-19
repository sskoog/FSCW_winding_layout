function kw = cros( Q,p,plot_on,verbose )
% Stefan Skoog 2017
% Chalmers University of Technology, Gothenburg, Sweden

% Function inputs:
% 'Q' = number of stator slot openings = number of teeth
% 'p' = number of rotor poles = number of (surface mount) magnets
% 'plot_on' = non-zero values will plot phasor and geometry graphics
% 'verbose' = non-zero values prints a lot of text in console
% Two example motors:
% Automotive starter-generator: p = 32, Q = 24
% RC motor for hobby airplanes: p = 20, Q = 18


% Implementation of Cros' method of finding the optimal winding layout for
% fractional-slot concentrated-winding PMSM.
% This algorithm will create the optimal winding layout for a two-layer
% FSCW and calculate the winding factor through phasor evaluation

% The methodology is explained in reserach papers:
% Source 1: Permanent-Magnet Synchronous Machines with Non-Overlapping
% Concentrated Windings for Low-Speed Direct-Drive Applications.
% Florence Meier, PhD Thesis 2008 KTH

% Source 2:
% Investigation on pole-slot combinations for permanent-magnet machines
% with concentrated windings"
% F Libert, J Soulard
% International Conference on Electrical Machines (ICEM2004)

% Source 3:
% Winding Factors and Joule Losses of Permanent Magnet Machines with Concentrated Windings
% Freddy Magnussen and Chandur Sadarangani

if nargin < 3
    plot_on = 1;
    verbose = 1;
elseif nargin < 4
    verbose = 1;
end

% Number of phases (only tested for 3)
n_ph = 3;

% Fallback return value
kw = 0;

% Reference winding for a 3-phase machine with full-pitch windings
% Capital letters are inflow, lower-case is return
ref_wind = 'AcBaCb';

if verbose fprintf('Cros'' method for fractional-slot pitch concentrated winding permanent magnetized synchronous machines\n'); end
if verbose fprintf('%i slots, %i poles, %i phases\n',Q,p,n_ph); end

assert(mod(p,2)==0,'Pole number must be an even number!');

%% Step A
% Calculate slot/pole/phase and reduce to fraction of lowest terms
% q = a/b
q = Q/p/n_ph;
[a b] = rat(q);
if verbose fprintf('%1.3f slots/pole/phase = %i/%i\n',q,a,b); end
if q>=1
	fprintf('Bad design @ Q/p=%i/%i: This method only works for fractional-slot machines: q<1 (q=%1.3f)!\n',Q,p,q);
    kw = NaN;
    return;
end

%% Step B
% Create sequence of 0s and 1s.
% a ones and |b-a| zeroes
% Distribute zeroes and ones as even as possible
ze = abs(b-a); % Number of zeros
seq = zeros(1,ze+a);
%ons = round(1:(ze+a)/a:(ze+a));
ons = floor(1:(ze+a)/a:(ze+a));
seq(ons) = 1;

if verbose
    fprintf('Sequence of %i ''0'' and %i ''1'':\n',ze,a);
    fprintf('    ');
    for i=1:length(seq)
        fprintf('%i ',seq(i));
    end
    fprintf('\n')
end

%% Step C
% Repeat sequence until all slots are populated
% That is: 3*p/b = Q/a
seq2 = seq;
for i=2:round(Q/a)
    seq2 = cat(2,seq2,seq);
end
if verbose
    fprintf('Repeating sequence %i times. Total length: %i. Compare with reference winding:\n',Q/a,length(seq2));
    fprintf('    ');
    for i=1:length(seq2)
        fprintf('%i ',seq2(i));
    end
    fprintf('\n')
end

%% Step D
% Compare sequence with the reference winding
% Keep reference winding layout only where '1' are present in sequence
% Each new '1' from the sequence is placed at a new tooth
% Then, insert return conductor at each tooth so that all slots are filled

% Create reference sequence with at leaste same length as sequence
ref2 = ref_wind;
while length(ref2) < length(seq2)
    ref2 = cat(2,ref2,ref_wind);
end
if verbose
    fprintf('    ');
    for i=1:length(seq2)
        fprintf('%s ',ref2(i));
    end
    fprintf('\n')
end
% Identify all '1's in sequence
idx = find(seq2==1);
% Pick all indeces of long reference vector that corresponds to '1's
wind1 = ref2(idx);

% if verbose
%     fprintf('Resulting half-sequence:\n  ');
%     for i=1:length(wind1)
%         fprintf('%s|*  ',wind1(i));
%     end
%     fprintf('\n')
% end

% Insert matching return winding
for i=1:length(wind1)
    c1 = wind1(i);
    c2 = returnWinding(wind1(i));
    wind2(2*i-1) = c1;
    wind2(2*i) = c2;
	%fprintf('  Input: %s. Output: %s|%s  (%i|%i)\n',c2,c1,c2,2*i-1,2*i);
end
if verbose
    fprintf('Full sequence and slot number:\n  ');
    i=1;
    while i<length(wind2)
        fprintf('%s|',wind2(i));
        i=i+1;
        fprintf('%s  ',wind2(i));
        i=i+1;
    end
    fprintf('\n')
    for i=1:round(length(wind2)/2)
        fprintf('%02.0f   ',i);
    end
    fprintf('\n')
end

%% Step E
% Calculate vector S, to describe the layout of phase A
% S is the base of winding factor calculation
% Name all slots in sequence from 1 to Q
% In S, map the index to all Phase A conductors. For return conductors,
% append a - to the index.
k1=1;
k2=1;
k3=1;
S1 = [];
S2 = [];
S3 = [];
for i=1:length(wind2)
    slot = floor(i/2)+1;
    % Phase A
    if wind2(i) == 'A'
        S1(k1) = slot;
        k1 = k1+1;
    elseif wind2(i) == 'a'
        S1(k1) = -slot;
        k1 = k1+1;
    end
    % Phase B
    if wind2(i) == 'B'
        S2(k2) = slot;
        k2 = k2+1;
    elseif wind2(i) == 'b'
        S2(k2) = -slot;
        k2 = k2+1;
    end
    % Phase C
    if wind2(i) == 'C'
        S3(k3) = slot;
        k3 = k3+1;
    elseif wind2(i) == 'c'
        S3(k3) = -slot;
        k3 = k3+1;
    end
end

% Design sanity check 1
if length(S1)<1 | length(S2)<1 | length(S3)<1
    fprintf('Bad design @ Q/p=%i/%i: No phase B or C windings!\n',Q,p);
    kw = NaN;
    return;
end

if verbose
    fprintf('S1 vector: ');
    for i=1:length(S1)
        fprintf('%i ',S1(i));
    end
    fprintf('\n')
    fprintf('S2 vector: ');
    for i=1:length(S2)
        fprintf('%i ',S2(i));
    end
    fprintf('\n')
    fprintf('S3 vector: ');
    for i=1:length(S3)
        fprintf('%i ',S3(i));
    end
    fprintf('\n')
end

% Design sanity check 2
if std([length(S1) length(S2) length(S3)]) > 0.1
    fprintf('Bad design @ Q/p=%i/%i: Imbalanced windings!\n',Q,p);
    kw = NaN;
    return;
end

%% Step F:
% Create Electro-Motive Force Phasors E from the S vector
%E = exp(S.*j*pi*p/Q);
E1 = sign(S1).*exp(abs(S1).*j*pi*p/Q);
E2 = sign(S2).*exp(abs(S2).*j*pi*p/Q);
E3 = sign(S3).*exp(abs(S3).*j*pi*p/Q);
% Calculate winding factor (same for all phases)
kw = abs(sum(E1))*3/(2*Q);
if verbose fprintf('Total winding factor: %1.3f \n',kw); end;

if plot_on
    % Pretty colormap. Creds to ColorBrewer 'set1'
    cm = {  [0.1059    0.6196    0.4667]...
            [0.8510    0.3725    0.0078]...
            [0.4588    0.4392    0.7020]};
    LineColorQ = [1 1 1]*0.4;
    LineColorP = [0.5 0.5 1]*0.6;

    f = figure();
    set(f,'Units','Normalized','OuterPosition',[0.0, 1-0.6, 0.30, 0.60]);
    % Create three axes
    % ax0 is for phasor illustration
    %ax0 = subplot(2,1,1);
    ax0 = axes('Position',[0.0 0.40 1 0.6]);
    set(ax0,'XTick',[]); % Disable x and y axes tiks
    set(ax0,'XTickLabel',[]);
    set(ax0,'ytick',[]);
    set(ax0,'yticklabel',[]);
    hold(ax0,'on');
    % ax1 is for stator tooth poles
    %ax1 = subplot(2,1,2);
    ax1 = axes('Position',[0.01 0.08 0.98 0.23]);
    ax1.XColor = LineColorQ;
    set(ax1,'ytick',[]);
    set(ax1,'yticklabel',[]);
    hold(ax1,'on');
    % Make another x-axis for magnet poles
    ax2 = axes('Position',ax1.Position,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    ax2.XColor = LineColorP;
    set(ax2,'ytick',[]);
    set(ax2,'yticklabel',[]);
    hold(ax2,'on');

    % Sweep phasors in ax0
    for i=1:length(E1)
        if i==1
            oldx1 = real(E1(i));
            oldy1 = imag(E1(i));
            quiver(ax0,0,0,oldx1,oldy1,0,'LineWidth',2); 
        else
            quiver(ax0,oldx1,oldy1,real(E1(i)),imag(E1(i)),0,'LineWidth',2);
            oldx1 = oldx1 + real(E1(i));
            oldy1 = oldy1 + imag(E1(i));
        end

        % Phase B
        if i==1
            oldx2 = real(E2(i));
            oldy2 = imag(E2(i));
            quiver(ax0,0,0,oldx2,oldy2,0,'LineWidth',2); hold on;
        else
            quiver(ax0,oldx2,oldy2,real(E2(i)),imag(E2(i)),0,'LineWidth',2);
            oldx2 = oldx2 + real(E2(i));
            oldy2 = oldy2 + imag(E2(i));
        end

        % Phase C
        if i==1
            oldx3 = real(E3(i));
            oldy3 = imag(E3(i));
            quiver(ax0,0,0,oldx3,oldy3,0,'LineWidth',2); hold on;
        else
            quiver(ax0,oldx3,oldy3,real(E3(i)),imag(E3(i)),0,'LineWidth',2);
            oldx3 = oldx3 + real(E3(i));
            oldy3 = oldy3 + imag(E3(i));
        end
    end
    quiver(ax0,0,0,oldx1,oldy1,0,'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])
    quiver(ax0,0,0,oldx2,oldy2,0,'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])
    quiver(ax0,0,0,oldx3,oldy3,0,'LineWidth',1.5,'LineStyle','--','Color',[0 0 0])
    axis(ax0,'equal');
    text(ax0,0.05*ax0.XLim(2),0.8*ax0.YLim(2),...
        sprintf('Q_s/p = %i/%i vector potentials',Q,p),...
        'FontSize',13,'FontWeight','bold');
    text(ax0,0.05*ax0.XLim(2),0.6*ax0.YLim(2),...
        sprintf('Winding factor = %1.3f',kw),...
        'FontSize',13,'FontWeight','bold');


    %% Step G: Draw slot/magnet geometry basic layout
    yp1 = 1;  % Start of magnet pole in y-dir
    yp2 = 5;  % End of magnet pole in y-dir
    yQ1 = -1; % Start of slot in y-dir
    yQ2 = -20; % End of slot in y-dir. Stator back.
    TP = 0.40; % ToothPitch as amount of slot pitch
    MP = 0.70; % MagnetPitch as amount of pole pitch

    % Just for the legend
    quiver(ax1, -1, 0,0,5,...
            'LineWidth',3,'Color',cm{1},'AutoScale','off','MaxHeadSize',0.25);
    quiver(ax1, -1, 0,0,5,...
            'LineWidth',3,'Color',cm{2},'AutoScale','off','MaxHeadSize',0.25)
    quiver(ax1, -1, 0,0,5,...
            'LineWidth',3,'Color',cm{3},'AutoScale','off','MaxHeadSize',0.25)    

    % Use slot number as primary x-grid. Pole numbers can be higher than slot number in some cases though!
    for x=1:Q
        % Draw line for one pole pitch at rotor back line
        line(x+[-1 +1]/2,[yQ2 yQ2],'LineWidth',10,...
                'Color',LineColorQ,'Parent',ax1); % Full horizontal line
        hold on;
        % Draw surface mount magnet as rectangle
        rectangle('Position',[x-TP/2 yQ2 TP yQ1-yQ2],...
                  'EdgeColor',LineColorQ,'FaceColor',LineColorQ,...
                  'Parent',ax1);
        % Add arrow to indicate magnetization direction with positive current
        % Add arrow to indicate magnetization direction with positive current
        dir = 1; % Positive direction
        phase = wind1(x)-'A'+1;
        if phase > 3
            phase = wind1(x)-'a'+1;
            dir = -1;
        end
        quiver(ax1, x, yQ1+(yQ2-yQ1)/2,0,0.5*(yQ1-yQ2)*dir,...
            'LineWidth',3,'Color',cm{phase},'AutoScale','off','MaxHeadSize',0.25);      

    end


    % Draw magnet poles. They have to align with the stator slots in total length!
    %pScale = Q/p;
    pScale = 1;
    for x=1:p
        % Draw line for one pole pitch at rotor back line
        line((x+[-1 +1]/2)*pScale,[yp2 yp2],'LineWidth',10,...
                'Color',LineColorP,'Parent',ax2); % Full horizontal line
        hold on;
        % Draw surface mount magnet as rectangle
        rectangle('Position',[(x-MP/2)*pScale yp1 (MP)*pScale yp2-yp1],...
                  'EdgeColor',LineColorP,'FaceColor',LineColorP,...
                  'Parent',ax2);

        % Add arrow to indicate magnetization direction with positive current
        quiver(ax2, x, yp1+(yp2-yp1)/2,0,0.5*(yp2-yp1)*sign(rem(x,2)-0.5),...
            'LineWidth',3,'Color','k','AutoScale','off','MaxHeadSize',10);      

    end

    % Make two x-scales with tick each integer
    x1 = 1:1:Q;
    x2 = 1:1:p;
    xlim(ax1,[0 Q]+0.5);
    xlim(ax2,[0 p]+0.5);
    ax1.XTick = x1;
    ax2.XTick = x2;
    % Adjust Y limitations to perfectly room entire layout
    ylim(ax1,[yQ2 yp2]*1.2);
    ylim(ax2,[yQ2 yp2]*1.2);
    grid(ax1,'on');
    xlabel(ax1,'Stator tooth number');
    xlabel(ax2,'Rotor pole number');
    legend(ax1,'Phase A','Phase B','Phase C','Location','South');
    end

end

%% Support functions
function y = returnWinding(x)
    switch x
        case 'A'
            y = 'a';
        case 'a'
            y = 'A';
        case 'B'
            y = 'b';
        case 'b'
            y = 'B';
        case 'C'
            y = 'c';
        case 'c'
            y = 'C';
        otherwise
            y = '*';
    end
end



% EOF