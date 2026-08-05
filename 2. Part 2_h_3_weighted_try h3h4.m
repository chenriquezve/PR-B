clear;
clc;
rng(1992)
%Entering values:

%Variable values (info. up to 23-07-26 (inclusive)):
p_beat= .5;
steps = 1000000; %Beats to simulate
capacity = 3; %Of simultaneous bay-service
beat = 10; %Unit of measurement minutes; 1 beat equals to X minutes
service_beats = 3;
date = "03-08-2026";

% Arrival representations
yes_arrival = [0 0 1];
%b737_arrival = [0 0 1];
no_arrival = [0 0 0];

%Fixed values, does not change:
hours_per_day = 24;
mins_per_hour = 60;

    %should i need "effective days of service"??

%Calculating from above:
p_beat_arrival = p_beat;
p_beat_non_arrival = 1 - (p_beat_arrival);
beats_per_day = (hours_per_day * mins_per_hour) / beat;
%---------------------------P1---------------------------------------------
fprintf("Last time info. updated: %s, where p_beat_: %.4f, and p_non_arrival= %.4f\n", date, p_beat, p_beat_non_arrival)

%Generating all possible states:
h = 3;
b = 3;
valid_states = [];
for p1 = 0:b
    for p2 = 0:b
        for p3 = 0:b
            candidate_state = [p1 p2 p3];
                if sum(candidate_state) <= capacity
                    valid_states(end+1, :) = candidate_state;
                end
            end
     end
end
n_states_raw = size(valid_states, 1);

fprintf("Number of possible states: %d\n", n_states_raw);
fprintf("Possible states generated OK.-\n")


%Trimming to reachable states (J(0:3,3)):
possible_states = valid_states; %56 raw states
reachable_states = zeros(1, h); %Starting from 'empty state'
state_to_check = 1;

while state_to_check <= size(reachable_states, 1)
    current_state = reachable_states(state_to_check, :);
    shifted_state = [current_state(2:end), 0]; %Shifting state
    possible_next_states = shifted_state; %no arrival
    candidate_state = shifted_state; %e190 arrival
    candidate_state(service_beats) = candidate_state(service_beats) + 1;

    if sum(candidate_state) <= capacity
        possible_next_states(end+1, :) = candidate_state;
    end
   
    for next = 1:size(possible_next_states, 1) %possible next states
        candidate_state = possible_next_states(next, :);
        is_possible = ismember(candidate_state, possible_states, "rows");
        reached_OK = ismember(candidate_state, reachable_states, "rows");
        if is_possible && ~reached_OK
            reachable_states(end+1, :) = candidate_state;
        end
    end
    state_to_check = state_to_check + 1; %move to the next state
end

valid_states = reachable_states;
n_states = size(valid_states, 1);

discarded_states = n_states_raw - n_states;

fprintf("Number of reachable states: %d. Hence, %d states were discarded.\n", n_states, discarded_states);

fprintf("Reachable states are:\n");
disp(valid_states);

fprintf("\n--------------------------P1 RUN OK--------------------------\n")
%---------------------------P2---------------------------------------------

P_trans_matrix = zeros(n_states, n_states);

%Initialising rejection vectors:
reject_arrival = zeros(n_states, 1);

for current_state_number = 1:n_states
    current_state = valid_states(current_state_number, :);
    shifted_state = [current_state(2:end), 0]; %Shift, then append 0

%---Event 1: No arrival
    next_state = shifted_state;
    next_state_number = find(all(valid_states == next_state, 2)); %Finding state
    if isempty(next_state_number)
        fprintf("'No arrival' state was not found:\n")
        disp(next_state)
    end

    P_trans_matrix(current_state_number, next_state_number) = (P_trans_matrix(current_state_number, next_state_number) + p_beat_non_arrival); %Adding "p_non_arrival" to Trans. matrix "P"

%---Event 2: an arrival?
    candidate_state = shifted_state;
    candidate_state(service_beats) = candidate_state(service_beats) + 1;
    if sum(candidate_state) <= capacity
        next_state = candidate_state; %Accepted
    else
        next_state = shifted_state; %Rejected
        reject_arrival(current_state_number) = p_beat;
    end
    next_state_number = find(all(valid_states == next_state, 2)); %Finding state
    if isempty(next_state_number)
        fprintf("Following 'e190' state was not found:\n")
        disp(next_state)
    end
    P_trans_matrix(current_state_number, next_state_number) = (P_trans_matrix (current_state_number, next_state_number) + p_beat); %Adding "p_beat_e190" to Trans. matrix "P"
end

fprintf("Transition matrix (P), with arrival probs. = Markov matrix generated OK.-\n")
%disp(P_trans_matrix)

fprintf("Row sums (should all be 1, otherwise check that row's probs):\n")
disp(sum(P_trans_matrix, 2).')


%Heatmap plot in white and red:
red_map = [ones(256,1), linspace(1,0,256)', linspace(1,0,256)'];
P_plot = [P_trans_matrix, nan(n_states,1);nan(1,n_states + 1)];
figure;
matrix_plot = pcolor(0:n_states, 0:n_states, P_plot);
matrix_plot.EdgeColor = [0.75 0.75 0.75];
matrix_plot.LineWidth = 0.3;
colormap(red_map);
clim([0 max(P_trans_matrix(:))]);
ax = gca;
ax.XAxisLocation = "top";
ax.YDir = "reverse";
ax.FontSize = 12;
ax.XTick = 0.5:1:n_states;
ax.YTick = 0.5:1:n_states;
state_labels = join(string(valid_states), " ", 2);
ax.XTickLabel = state_labels;
ax.YTickLabel = state_labels;
axis tight;
axis square;
xlabel("Following state", "FontSize", 14);
ylabel("Current state","FontSize", 14);
title("Transition probability matrix 'P' for J(0:3,3)","FontSize", 20);
c = colorbar;
c.FontSize = 14;
c.Label.String = "Transition probability";
xtickangle(ax,0);
ytickangle(ax,0);
ax.XTickLabelRotationMode = "manual";
fprintf("\n--------------------------P2 RUN OK--------------------------\n")
%---------------------------P3---------------------------------------------

%Step 3 - Checking whether a distribution (vector v) is stationary for P (transition matrix):
    [eig_vector, eig_value] = eig (P_trans_matrix');
    eigenvalues = diag(eig_value); %extracting eigenvalue
    fprintf("Eigenvalues:\n")
    disp(diag(eig_value'))
    tolerance = 0.00001;
    column_eig = find(abs(eigenvalues - 1) < tolerance); %Automatically choosing the column with the value closer to 1
    fprintf("Number of eigenvalues equal to 1: %d\n", length(column_eig))
    if length(column_eig) ~= 1
        fprintf("There is %d egenvalues equal to 1!", length(column_eig))
    end
    fprintf("Selected eigenvector column: %d\n", column_eig)
    v = eig_vector(:, column_eig); % extracting eigenvector
    v = v / sum(v); %Normalising the eigenvector (so the probs add to 1)
    fprintf("Double checking (v * P = v)") %Double checking selection is correct
    v_step = v' * P_trans_matrix;
    diff_v = norm(v_step - v'); % Difference between v and v*P (should be close to zero)
    fprintf("    Difference (should be close to 0): %.10f\n", diff_v)
    fprintf("    Used tolerance: %.4d\n", tolerance);
    if diff_v < tolerance
        fprintf("    -> v IS the stationary distribution\n")
    else
        fprintf("    -> v is NOT the stationary distribution, please check it again\n")
    end

    fprintf("Theoretical stationary distribution:\n")
    disp(v')
    

    number_eigenvalues = length(eigenvalues);
    eigenvalue_modulo = abs(eigenvalues);
    figure;
    eigenvalue_bars = bar(1:number_eigenvalues, eigenvalue_modulo);
    xlabel("Eigenvector column", "FontSize", 14);
    ylabel("Absolute eigenvalue (|\lambda|)", "FontSize", 14);
    title("J(0:3,3) - Absolute values of the transition matrix eigenvalues", "FontSize", 16);
    xticks(1:number_eigenvalues);
    ylim([0, 1]);

fprintf("\n--------------------------P3 RUN OK--------------------------\n")
%---------------------------P4---------------------------------------------
%Iteration to check convergence to the stationary distribution:
%Starting from a single state, randomly jumping to neighbours
% Recording visit frequency to approximate stationary distribution


%Plus, Im going to plot to visually represent the distribution convergence,
%the distance will be measured using residual sum of squares (RSS), which
%shows how far the current position (or current distribution) is from the
%theoretical stationary ditribution (v).
%Intuitively the smaller RSS, the closer to v.

fprintf("For h=3; capacity=3:\n")
%---------------------Starting point: p0(1)--------------------------------
fprintf("\n________________________SIMULATION 1________________________\n")
fprintf("Starting from State (p0(1)):\n") 
curr_state = find(all(valid_states == [0 0 0], 2), 1); %Starting from [0 0 0]
visit_count = zeros (1, n_states);%To store state visits at each step
rss_time_state1 = zeros (steps, 1); %To store RSS at each step
state_history1 = zeros(steps,1); %To store States at each step
simulation_counter1 = 0; %Seated e190
rejected_counter1 = 0; %Rejected 737
fprintf("    Set starting state:")
disp(curr_state)
fprintf("    Visit counts set to 0.-")

for time = 1:steps
    prev_state = curr_state;
    current_vec = valid_states(prev_state, :);
    shifted_vec = [current_vec(2:end), 0]; %Depart + shift
    r = rand;
    if r < p_beat_non_arrival %Event 1: No arrival
        next_vec = shifted_vec;
    elseif r < p_beat_non_arrival + p_beat% Event 2: arrival
        candidate_vec = shifted_vec;
        candidate_vec(service_beats) = candidate_vec(service_beats) + 1;
        if sum(candidate_vec) <= capacity
            next_vec = candidate_vec; %Seated
            simulation_counter1 = simulation_counter1 + 1;
        else
            next_vec = shifted_vec; %Rejected
            rejected_counter1 = rejected_counter1 + 1;
        end
    end
    
    curr_state = find(all(valid_states == next_vec, 2), 1);
    state_history1(time) = curr_state;
    visit_count(curr_state) = visit_count(curr_state) +1;
    observed_distr = visit_count/time;
    rss_time_state1(time) = sum ((observed_distr - v').^2);
end

sim_seated_total1 = simulation_counter1;
sim_rejected_total1 = rejected_counter1;
sim_arrivals_total1 = sim_seated_total1 + sim_rejected_total1;
sim_rejection_rate1 = (rejected_counter1/(simulation_counter1+rejected_counter1))*100;
sim_rejection_rate_total1 = (sim_rejected_total1/sim_arrivals_total1)*100;
sim_arrivals1 = simulation_counter1+ rejected_counter1;

sim_arr_beat1  = sim_arrivals1/steps;
sim_serv_beat1= simulation_counter1/steps;
sim_rej_beat1 = rejected_counter1/steps;


fprintf("\n\nRESULTS - From State 1 [0 0 0] (empty):\n")

fprintf("After %d beats:\n", steps)
fprintf(" -e190:\n")
fprintf("    Arrivals: %d\n", sim_arrivals1)
fprintf("    Effectively seated: %d\n", simulation_counter1)
fprintf("    Rejected: %d\n", rejected_counter1)


fprintf("\nPer beat:\n")
fprintf(" \n    Arrivals %.4f\n    Effectively seated %.4f\n    Rejected: %.4f\n", sim_arr_beat1, sim_serv_beat1, sim_rej_beat1)


fprintf("\n    Rejection rate total: %.4f%%\n", sim_rejection_rate_total1)
fprintf("    Rejection rate e190: %.4f%%\n", sim_rejection_rate1)

%---------------------Starting point: p0(8)--------------------------------
fprintf("\n")
fprintf("\n________________________SIMULATION 2________________________\n")
%	State 8	->	[1 1 1] 

curr_state = find(all(valid_states == [1 1 1], 2), 1); %Starting from [1 1 1]
fprintf("Starting from State (p0(%d)):\n", curr_state)
visit_count = zeros (1, n_states); %To count states at each visit
rss_time_state2 = zeros (steps, 1); %To store RSS at each step
state_history2 = zeros(steps,1); %To store States at each step
simulation_counter2 = 0; %Seated
rejected_counter2 = 0; %Rejected

fprintf("    Set starting state:")
disp(curr_state)
fprintf("    Visit counts set to 0.-")

for time = 1:steps
prev_state = curr_state;
    current_vec = valid_states(prev_state, :);
    shifted_vec = [current_vec(2:end), 0];   %Depart + shift
    r = rand;
    if r < p_beat_non_arrival %Event 1: No arrival
        next_vec = shifted_vec;
    elseif r < p_beat_non_arrival + p_beat%Event 2:arrival
        candidate_vec = shifted_vec;
        candidate_vec(service_beats) = candidate_vec(service_beats) + 1;
        if sum(candidate_vec) <= capacity
            next_vec = candidate_vec; %Seated
            simulation_counter2 = simulation_counter2 + 1;
        else
            next_vec = shifted_vec; %Rejected
            rejected_counter2 = rejected_counter2 + 1;
        end
    end
 
    curr_state = find(all(valid_states == next_vec, 2), 1);
    state_history2(time) = curr_state;
    visit_count(curr_state) = visit_count(curr_state) + 1;
    observed_distr = visit_count/time;
    rss_time_state2(time) = sum((observed_distr - v').^2);
end

sim_seated_total2 = simulation_counter2;
sim_rejected_total2 = rejected_counter2;
sim_arrivals_total2 = sim_seated_total2 + sim_rejected_total2;
sim_rejection_rate2 = (rejected_counter2/(simulation_counter2+rejected_counter2))*100;
sim_rejection_rate_total2 = (sim_rejected_total2/sim_arrivals_total2)*100;
sim_arrivals2 = simulation_counter2 + rejected_counter2;

sim_arr_beat2 = sim_arrivals2/steps;
sim_serv_beat2 = simulation_counter2/steps;
sim_rej_beat2 = rejected_counter2/steps;

sim_arr_beat2_total  = sim_arrivals_total2/steps;
sim_serv_beat2_total = sim_seated_total2/steps;
sim_rej_beat2_total  = sim_rejected_total2/steps;

fprintf("\n\nRESULTS - From State 8 [1 1 1] (full):\n")

fprintf("After %d beats:\n", steps)
fprintf(" :\n")
fprintf("    Arrivals: %d\n", sim_arrivals2)
fprintf("    Effectively seated: %d\n", simulation_counter2)
fprintf("    Rejected: %d\n", rejected_counter2)

fprintf(" -Total:\n")
fprintf("    Arrivals: %d\n", sim_arrivals_total2)
fprintf("    Serviced: %d\n", sim_seated_total2)
fprintf("    Rejected: %d\n", sim_rejected_total2)


fprintf("\nPer beat:\n")
fprintf(" \n    Arrivals %.4f\n    Effectively seated %.4f\n    Rejected: %.4f\n", sim_arr_beat2, sim_serv_beat2, sim_rej_beat2)
fprintf(" -Total:\n    Arrivals %.4f\n    Effectively seated %.4f\n    Rejected: %.4f\n", sim_arr_beat2_total, sim_serv_beat2_total, sim_rej_beat2_total)


fprintf("\n    Rejection rate total: %.4f%%\n", sim_rejection_rate_total2)
fprintf("    Rejection rate: %.4f%%\n", sim_rejection_rate2)

fprintf("\n___________RSS DIFF. VALUES FOR MY REPORT___________\n")
checkpoints = [10, 1000, 100000, steps];
fprintf("Step\t\tEmpty start\t\tFull start\n");
for checkpoint = 1:length(checkpoints)
    fprintf("%d\t\t%.6e\t\t%.6e\n", checkpoints(checkpoint), rss_time_state1(checkpoints(checkpoint)), rss_time_state2(checkpoints(checkpoint)));
end



%------------Warm up time together in log: p0(1) & p0(39)------------------

figure;
loglog(1:steps, rss_time_state1, 'LineWidth', 1.2); % empty start
hold on;
loglog(1:steps, rss_time_state2, 'LineWidth', 1.2); % full start
hold off;
xlabel("Number of steps - log scale", "FontSize", 14);
ylabel("Residual Sum of Squares (RSS) - log scale", "FontSize", 14);
title("J(0:3,3) - Warm up time comparison: empty vs. full start", "FontSize", 16);
legend("Started empty [0 0 0]", "Started full [1 1 1]", "Location", "northeast", "FontSize", 12);

%-----------------Both together - visit freq: p0(1) & p0(8)---------------

freq_from_empty = zeros(n_states, 1); %starting from [0 0 0]
freq_from_full = zeros(n_states, 1); %starting from [1 1 1]

%Normalising:
for state = 1:n_states
    freq_from_empty(state) = sum(state_history1 == state) / steps;
    freq_from_full(state)  = sum(state_history2 == state) / steps;
end

state_labels = join(string(valid_states), "", 2);
comparison_values = [v, freq_from_empty, freq_from_full];

[ranked_v, ranked_order] = sort(v, "descend");
ranked_state_labels = state_labels(ranked_order);
ranked_freq_empty = freq_from_empty(ranked_order);
ranked_freq_full = freq_from_full(ranked_order);
ranked_comparison_values = [ranked_v, ranked_freq_empty, ranked_freq_full];
cumulative_v = cumsum(ranked_v) * 100;
cumulative_empty = cumsum(ranked_freq_empty) * 100;
cumulative_full = cumsum(ranked_freq_full) * 100;

%Ranking all states
figure;
bar(ranked_comparison_values);
xticks(1:n_states);
xticklabels(ranked_state_labels);
xtickangle(90);
xlabel("States ranked by theoretical probability", "FontSize", 14);
ylabel("Long-run probability", "FontSize", 14);
title("J(0:3,3) - Theoretical vs. simulated frequencies ", "FontSize", 16);
legend("Theory (v)", "Simulated - started empty [0 0 0]", "Simulated - started full [1 1 1]", "Location", "northeast", "FontSize", 12);

%Ranked values:
state_labels = join(string(valid_states), "", 2);
fprintf("\nRANKING OF STATES\n")
fprintf("Rank   State   Theory   Empty start    Full start    Cum. theory     Cum. empty     Cum. full     \n")

for rank_number = 1:min(10, n_states)
    fprintf("%d\t%s\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n", ...
        rank_number, char(ranked_state_labels(rank_number)), ranked_v(rank_number) * 100,...
        ranked_freq_empty(rank_number) * 100, ranked_freq_full(rank_number) * 100, ...
        cumulative_v(rank_number), cumulative_empty(rank_number), ...
        cumulative_full(rank_number));
end

fprintf("\n--------------------------P4 RUN OK--------------------------\n")


%---------------------------P5---------------------------------------------
%========== (1) Theoretical occupancy =========================
%Checking theoretical long-term bay occupancy to answer: "how full will the system be?":
p_empty = 0; %starting counter of 0 bays busy
p_busy_bay_1 = 0; %starting counter of 1 bay busy
p_busy_bay_2 = 0; %starting counter of 2 bays busy
p_full = 0; %starting counter of 3 bays busy

for state = 1:n_states
    occupancy_now = sum(valid_states(state, :));
    if occupancy_now == 0
        p_empty = p_empty + v(state);
    elseif occupancy_now == 1
        p_busy_bay_1 = p_busy_bay_1 + v(state);
    elseif occupancy_now == 2
        p_busy_bay_2 = p_busy_bay_2 + v(state);
    elseif occupancy_now == 3
        p_full = p_full + v(state);
    end
end

double_check_prob = p_empty + p_busy_bay_1 + p_busy_bay_2 + p_full;

fprintf("Theoretical long-run occupancy shares (from v):\n")
fprintf("   P(all bays empty)  = %.4f\n", p_empty)
fprintf("   P(1 bay busy)     = %.4f\n", p_busy_bay_1)
fprintf("   P(2 bays busy)    = %.4f\n", p_busy_bay_2)
fprintf("   P(full)   = %.4f\n", p_full)
fprintf("   Check (should sum to 1): %.4f\n", double_check_prob)

mean_busy_bays = 0*p_empty + 1*p_busy_bay_1 + 2*p_busy_bay_2 + 3*p_full;
utilisation = mean_busy_bays / capacity;
fprintf("   Mean number of busy bays = %.4f of %d\n", mean_busy_bays, capacity)
fprintf("   Long-run utilisation = %.2f%%\n", utilisation * 100)

%========== (2) Theoretical serviced and rejected =========================
%Counting planes that are being, theoretically in the long term, served and rejected:

served_per_beat = 0; %starting counter from 0
rejected_per_beat = 0; %starting counter from 0

for state = 1:n_states
    shifted_state = [valid_states(state, 2:end), 0];
    candidate_state = shifted_state;
    candidate_state(service_beats) = candidate_state(service_beats) + 1;
    if sum(candidate_state) <= capacity
        served_per_beat = served_per_beat + v(state) * p_beat;
    else
        rejected_per_beat = rejected_per_beat + v(state) * p_beat;
    end
end

%Per beat:
%---Both:
served_total_per_beat   = served_per_beat;
rejected_total_per_beat = rejected_per_beat;

%Per day counters:
%---e190:
serviced_per_day = served_per_beat * beats_per_day;
rejected_per_day = rejected_per_beat * beats_per_day;
arrival_per_day = p_beat* beats_per_day;

fprintf("\nThe expected parameters in the theoretical long run p/plane-p/beat are as follows:\n")
fprintf("Note: arrival rates were estimated from ADL Airport data using Python.\n")
fprintf("Aircraft\tArrivals\tRejected\tServiced\n")
fprintf(" \t\t%.4f\t\t%.4f\t\t%.4f\n", p_beat, rejected_per_beat, served_per_beat);
fprintf("Total\t\t%.4f\t\t%.4f\t\t%.4f\n", p_beat_arrival, rejected_total_per_beat, served_total_per_beat);


fprintf("\n\nThe expected parameters in the theoretical long run p/plane-p/day are as follows:\n")
fprintf("Aircraft\tArrivals\tRejected\tServiced\n")

fprintf("E190\t\t%.4f\t\t%.4f\t\t%.4f\n", arrival_per_day, rejected_per_day, serviced_per_day);

fprintf("\n--------------------------P5 RUN OK--------------------------\n")


%---------------------------P6--------------------------------------------
%========== (1) Occupancy: theory vs, simulated =========================
%Theory (p_empty..p_full) from P6. Simulated summed from freq_from_empty.
sim_occ = zeros(4, 1); %positions = 0,1,2,3 bays busy
for state = 1:n_states
    occupancy_now = sum(valid_states(state,:));
    sim_occ(occupancy_now + 1) = sim_occ(occupancy_now + 1) + freq_from_empty(state);
end
theo_occ = [p_empty; p_busy_bay_1; p_busy_bay_2; p_full];

figure;
bar(0:3,[theo_occ, sim_occ]); %two bars per occupancy level
xlabel("Number of busy bays", "FontSize", 14);
ylabel("Long-term probability", "FontSize", 14);
title("J(0:3,3) - System occupancy levels: theory (v) vs. simulated (started empty)", "FontSize", 16);
legend("Theory (v)", "Simulated - started empty [0 0 0]", "FontSize", 12);

%Puede ser un dato por tabla

fprintf("Plot generated OK.-\n--------------------------P6 RUN OK--------------------------\n")

%---------------------------P7--------------------------------------------
%Sensitivity: how does rejection rate change if arrivals increasing?

next_N = zeros(n_states, 1);%after "no arrival"
next_S = zeros(n_states, 1);%after an arrival
for state = 1:n_states
    current_vec = valid_states(state, :);
    shifted_vec = [current_vec(2:end), 0];
    next_N(state) = find(all(valid_states == shifted_vec, 2), 1);%no arrival
    candidate_vec = shifted_vec;%arrival
    candidate_vec(service_beats) = candidate_vec(service_beats) + 1;
    if sum(candidate_vec) > capacity
        candidate_vec = shifted_vec;%rejected
    end
    next_S(state) = find(all(valid_states == candidate_vec, 2), 1);
end

% Arrival-probability amplification
factor_min = 0;
factor_step = 0.1;
factor_max = 2;
amplified = factor_min:factor_step:factor_max;
x_ampli = zeros(length(amplified), 1);
y_ampli = zeros(length(amplified), 1);

for m = 1:length(amplified)
    p_amplified = p_beat * amplified(m);
    p_n_amplified = 1 - p_amplified;
    P_trans_matrix_amplified = zeros(n_states, n_states);
    for state = 1:n_states
        P_trans_matrix_amplified(state, next_N(state)) = P_trans_matrix_amplified(state, next_N(state)) + p_n_amplified;
        P_trans_matrix_amplified(state, next_S(state)) =P_trans_matrix_amplified(state, next_S(state))+ p_amplified;
    end
    % Stationary distribution for the amplified model
    [eig_vec_ampli, eig_val_ampli] = eig(P_trans_matrix_amplified');
    eigenvalues_ampli = diag(eig_val_ampli);
    [~, col] = min(abs(eigenvalues_ampli - 1));
    v_ampli = real(eig_vec_ampli(:, col));
    v_ampli = v_ampli / sum(v_ampli);

    % Long-run probability that an arrival is rejected
    reject_prob = 0;
    for state = 1:n_states
        shifted_vec = [valid_states(state, 2:end), 0];
        candidate_vec = shifted_vec;
        candidate_vec(service_beats) = candidate_vec(service_beats) + 1;
        if sum(candidate_vec) > capacity
            reject_prob = reject_prob + v_ampli(state);
        end
    end
    x_ampli(m) = p_amplified;
    y_ampli(m) = reject_prob * 100;
end

valid = ~isnan(y_ampli);
figure;
plot(amplified(valid), y_ampli(valid), "-o","LineWidth", 1.5, "MarkerSize", 5);
hold on;
xline(1, "--", "Original estimated arrival probability", "LabelVerticalAlignment", "bottom", "FontSize", 12);
hold off;
xlabel("Arrival-probability multiplier", "FontSize", 14);
ylabel("Theoretical rejection percentage (%)", "FontSize", 14);
title("J(0:3,3) - Effect of amplified arrivals on rejection rate", "FontSize", 16);


fprintf("\n--------------------------P9 RUN OK--------------------------\n")

save("h3_remake.mat", "h", "valid_states", "n_states", ...
    "P_trans_matrix", "eig_vector", "eig_value", "v", "column_eig", ...
    "v_step", "diff_v")

disp("Correctly saved")
