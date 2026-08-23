% Exercise 1: Random DNA sequence generation and base frequency
rng(1);
bases = 'ATGC';
seq_len = 10000;
idx = randi(4, 1, seq_len);
dna = bases(idx);

counts = zeros(1, 4);
for i = 1:4
    counts(i) = sum(dna == bases(i));
end

figure;
bar(counts, 'FaceColor', [0.3 0.5 0.8]);
set(gca, 'XTickLabel', {'A', 'T', 'G', 'C'});
xlabel('Base');
ylabel('Count');
title('Base Frequency in Random DNA Sequence');
grid on;

fprintf('A: %d (%.2f%%)\n', counts(1), counts(1)/seq_len*100);
fprintf('T: %d (%.2f%%)\n', counts(2), counts(2)/seq_len*100);
fprintf('G: %d (%.2f%%)\n', counts(3), counts(3)/seq_len*100);
fprintf('C: %d (%.2f%%)\n', counts(4), counts(4)/seq_len*100);






% Exercise 2: GC content in sliding windows
window_size = 100;
n_windows = seq_len - window_size + 1;
gc_content = zeros(1, n_windows);

is_gc = (dna == 'G') | (dna == 'C');

for i = 1:n_windows
    window = is_gc(i:i+window_size-1);
    gc_content(i) = sum(window) / window_size;
end

figure;
plot(gc_content, 'LineWidth', 1);
yline(0.5, 'r--', 'LineWidth', 1.5, 'Label', 'Expected (50%)');
xlabel('Position (window start)');
ylabel('GC Content');
title('GC Content Along the Sequence (100-base sliding window)');
grid on;






% Exercise 3: Motif search
motif = 'ATAT';
motif_len = length(motif);
count = 0;

for i = 1:(seq_len - motif_len + 1)
    if strcmp(dna(i:i+motif_len-1), motif)
        count = count + 1;
    end
end

expected = seq_len * (0.25)^motif_len;

fprintf('Motif "%s" found %d times\n', motif, count);
fprintf('Expected count for random DNA: %.2f\n', expected);
fprintf('Ratio (observed/expected): %.2f\n', count/expected);







% Exercise 4: Match percentage -- random vs mutated
rng(2);
n = 1000;
idx1 = randi(4, 1, n);
idx2 = randi(4, 1, n);
seq1 = bases(idx1);
seq2 = bases(idx2);

match_random = sum(seq1 == seq2) / n * 100;
fprintf('Match %% between two random sequences: %.2f%%\n', match_random);

% Simulate evolution: mutate seq1
mutation_rates = [0.001, 0.01, 0.05, 0.1, 0.2, 0.5];
match_percentages = zeros(size(mutation_rates));

for k = 1:length(mutation_rates)
    mut_rate = mutation_rates(k);
    mutated = seq1;
    for i = 1:n
        if rand() < mut_rate
            other_bases = bases(bases ~= seq1(i));
            mutated(i) = other_bases(randi(3));
        end
    end
    match_percentages(k) = sum(seq1 == mutated) / n * 100;
end

figure;
plot(mutation_rates*100, match_percentages, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Mutation Rate (%)');
ylabel('Match Percentage (%)');
title('Mutation Rate vs Sequence Match Percentage');
grid on;

fprintf('\nMutation rate vs Match %%:\n');
for k = 1:length(mutation_rates)
    fprintf('  %.1f%% mutations -> %.2f%% match\n', mutation_rates(k)*100, match_percentages(k));
end







% Exercise 5: Markov chain DNA generator
% Using a sample "real" sequence as a stand-in for downloaded genomic data
% (replace this with fastaread('yourfile.fasta') if you have real data)
real_seq = ['ATGCGATCGTAGCTAGCTAGGGCATCGATCGATCGTAGCTAGCATCGATCGTAGCATGCATCG' ...
    'GATCGATCGGGGCATCGATCGCGCGCGATCGTAGCTAGCTAGCATGCATCGATCGGGGCATC' ...
    'GATCGATCGTAGCATGCATCGGGGCATCGATCGATCGTAGCTAGCATGCATCGATCGGGGCA'];

bases = 'ATGC';
n_bases = 4;
transition_counts = zeros(n_bases, n_bases);

for i = 1:length(real_seq)-1
    curr = strfind(bases, real_seq(i));
    next = strfind(bases, real_seq(i+1));
    transition_counts(curr, next) = transition_counts(curr, next) + 1;
end

transition_matrix = transition_counts ./ sum(transition_counts, 2);

disp('Transition matrix P(next base | current base):');
disp(array2table(transition_matrix, 'VariableNames', {'A','T','G','C'}, ...
    'RowNames', {'A','T','G','C'}));

% Generate artificial DNA using the transition matrix
n_generate = 1000;
markov_seq = zeros(1, n_generate);
markov_seq(1) = randi(4);  % random start

for i = 2:n_generate
    probs = transition_matrix(markov_seq(i-1), :);
    markov_seq(i) = find(rand() <= cumsum(probs), 1);
end

markov_dna = bases(markov_seq);

% Compare GC content: pure random vs Markov-generated vs real
gc_real = sum(real_seq == 'G' | real_seq == 'C') / length(real_seq) * 100;
gc_markov = sum(markov_dna == 'G' | markov_dna == 'C') / length(markov_dna) * 100;
gc_random = sum(dna == 'G' | dna == 'C') / length(dna) * 100;

fprintf('\nGC content comparison:\n');
fprintf('Real sequence:      %.2f%%\n', gc_real);
fprintf('Markov-generated:   %.2f%%\n', gc_markov);
fprintf('Pure random:        %.2f%%\n', gc_random);











% Exercise 6: Phylogenetic distance from a common ancestor
rng(3);
ancestor_len = 500;
ancestor = bases(randi(4, 1, ancestor_len));

% Create 4 "species" by mutating the ancestor at different rates
% (simulating different amounts of evolutionary time / divergence)
mutation_rates = [0.02, 0.05, 0.08, 0.15];
species_names = {'SpeciesA', 'SpeciesB', 'SpeciesC', 'SpeciesD'};
species_seqs = cell(1, 4);

for s = 1:4
    seq = ancestor;
    for i = 1:ancestor_len
        if rand() < mutation_rates(s)
            other_bases = bases(bases ~= ancestor(i));
            seq(i) = other_bases(randi(3));
        end
    end
    species_seqs{s} = seq;
end

% Compute pairwise distance matrix (% different bases)
n_species = 4;
dist_matrix = zeros(n_species, n_species);

for i = 1:n_species
    for j = 1:n_species
        diffs = sum(species_seqs{i} ~= species_seqs{j});
        dist_matrix(i, j) = diffs / ancestor_len * 100;
    end
end

disp('Pairwise distance matrix (% different bases):');
disp(array2table(dist_matrix, 'VariableNames', species_names, 'RowNames', species_names));

figure;
imagesc(dist_matrix);
colorbar;
colormap('hot');
set(gca, 'XTick', 1:4, 'XTickLabel', species_names, 'YTick', 1:4, 'YTickLabel', species_names);
title('Phylogenetic Distance Matrix');

% Find most closely related pair (excluding diagonal)
dist_matrix_nodiag = dist_matrix + diag(inf(1, n_species));
[min_val, min_idx] = min(dist_matrix_nodiag(:));
[row, col] = ind2sub(size(dist_matrix_nodiag), min_idx);
fprintf('\nMost closely related: %s and %s (%.2f%% different)\n', ...
    species_names{row}, species_names{col}, min_val);












