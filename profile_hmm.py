import argparse
import json
import numpy as np
import read_alignment as ra


TRANSITIONS = ("BM", "BI", "BD", "MM", "MI", "MD", "IM",
               "II", "ID", "DM", "DI", "DD", "ME", "IE", "DE")

BEGIN_TRANSITIONS = ("BM", "BI", "BD", "IM", "II", "ID")
MIDDLE_TRANSITIONS = ("MM", "MI", "MD", "IM", "II", "ID", "DM", "DI", "DD")
END_TRANSITIONS = ("MI", "ME", "II", "IE", "DI", "DE")

BACKGROUND_FREQUENCIES = "aa_frequencies.json"

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class ProfileHMM:
    def __init__(self, multiple_alignment, pseudocount=0.1):
        self.num_match_states, _ = get_match_states(multiple_alignment)
        self.transition_matrix = get_transition_matrix(
            multiple_alignment, pseudocount)
        self.m_emission_matrix = get_emission_matrix(
            multiple_alignment, pseudocount)
        self.i_emission_matrix = {}

        with open(BACKGROUND_FREQUENCIES, 'r') as file:
            self.i_emission_matrix = json.load(file)

    def align_sequence(self, sequence):
        """
        Given a sequence, returns the most probable path through the profile HMM.

        Args:
            - sequence: a string representing a sequence

        Returns:
            - path: a string representing the most probable path through the profile HMM
        """
        path = ""
        return path


def get_alignment_columns(multiple_alignment):
    """
    Given a multiple alignment, returns a list of columns, where each column is a list of characters, one for each sequence in the alignment

    Args:
        - multiple_alignment: a list of strings representing the sequences in the multiple alignment

    Returns:
        - columns: a list of columns, where each column is a string of characters, one for each sequence in the alignment
    """
    columns = []
    for i in range(len(multiple_alignment[0])):
        column = ""
        for sequence in multiple_alignment:
            column += sequence[i]
        columns.append(column)
    return columns


def get_alignment_columns_with_match_states(multiple_alignment):
    """
    Same as get_alignment_columns, but replaces dashes in non-match alignment columns with pound signs. This is to distinguish gaps due to deletion (in match states) and gaps due to insertion (in non-match states).

    Args:
        - multiple_alignment: a list of strings representing the sequences in the multiple alignment

    Returns:
        - columns: a list of columns, where each column is a string of characters, one for each sequence in the alignment. Dashes in non-match columns are replaced with pound signs.

    Example:
        >>> get_alignment_columns_with_match_states(['ABCE', 'AB-D', 'A--E', 'AC-D'])
        ['AAAA', 'BB-C', 'C###', 'EDED']
    """
    match_states, match_pattern = get_match_states(multiple_alignment)
    alignment_columns = get_alignment_columns(multiple_alignment)

    new_alignment_columns = []
    for col, match in zip(alignment_columns, match_pattern):
        if match:
            new_alignment_columns.append(col)
        else:
            new_alignment_columns.append(col.replace('-', '#'))

    return new_alignment_columns


def get_match_states(multiple_alignment):
    """
    Given a multiple alignment, returns the number of match states in the profile HMM and a string with the same length as the multiple alignment, with a * for match states and a space for other states

    Args:
        - multiple_alignment: a list of strings representing the sequences in   the multiple alignment

    Returns:
        - match_states: the number of match states in the profile HMM
        - match_pattern: a list of booleans, where True indicates a match state and False indicates a non-match state
    """
    num_sequences = len(multiple_alignment)
    match_states = 0
    match_pattern = []
    for col in range(len(multiple_alignment[0])):
        count_gaps = 0
        for seqeuence in multiple_alignment:
            if seqeuence[col] == "-":
                count_gaps += 1
        if count_gaps < num_sequences // 2:
            match_states += 1
            match_pattern.append(True)
        else:
            match_pattern.append(False)
    return (match_states, match_pattern)


def get_transition_counts(multiple_alignment):
    """
    Given a multiple alignment, returns a list of dictionaries, where each dictionary contains the number of times each type of transition appears each column corresponding to a match state.
    """
    num_sequences = len(multiple_alignment)
    match_states, match_pattern = get_match_states(multiple_alignment)
    alignment_columns = get_alignment_columns_with_match_states(
        multiple_alignment)

    # Group alignment columns by match states
    match_groups = [[] for _ in range(match_states)]
    match_index = -1
    for col, match in enumerate(match_pattern):
        if match:
            match_index += 1
            if col > 0:
                match_groups[match_index-1].append(alignment_columns[col])
        match_groups[match_index].append(alignment_columns[col])

    # Add begin and end states
    match_groups = [["^" * num_sequences, match_groups[0][0]]] + match_groups
    match_groups[-1].append("$" * num_sequences)

    # Separate match groups into transitions
    transitions = []
    for group in match_groups:
        transitions.append([''.join(e) for e in list(zip(*group))])

    transition_counts = [
        get_transition_counts_helper(t) for t in transitions]

    return transition_counts


def get_transition_counts_helper(transitions):
    """
    Given a list of transitions, returns a dictionary with the number of times each transition appears in the list.

    Args:
        - transitions: a list containing strings whereby each pair of characters represents a transition

    Returns:
        - transition_counts: a dictionary with the number of times each type of transition appears in the list
    """
    transition_counts = {t: 0 for t in TRANSITIONS}
    for string in transitions:
        clean_str = string.replace("#", "")
        str_len = len(clean_str)
        if str_len == 2:
            key = ('B' if clean_str[0] == '^' else 'D' if clean_str[0] == '-' else 'M') + \
                ('E' if clean_str[1] ==
                 '$' else 'D' if clean_str[1] == '-' else 'M')
            transition_counts[key] += 1
        else:
            pairs = [clean_str[i:i+2] for i in range(len(string) - 1)]

            start_pair = pairs.pop(0)
            start_key = (
                'B' if start_pair[0] == '^' else 'D' if start_pair[0] == '-' else 'M') + 'I'
            transition_counts[start_key] += 1

            end_pair = pairs.pop(-1)
            end_key = 'I' + \
                ('E' if end_pair[1] ==
                 '$' else 'D' if end_pair[1] == '-' else 'M')
            transition_counts[end_key] += 1

            for _ in pairs:
                transition_counts['II'] += 1
    return transition_counts


def get_transition_matrix(multiple_alignment, pseudocount=0.1):
    """
    Given a multiple alignment, returns a list of dictionaries, where each dictionary contains the transition probabilities for each column corresponding to a match state.
    """
    transition_counts = get_transition_counts(multiple_alignment)
    match_states, match_pattern = get_match_states(multiple_alignment)

    transition_matrix = [{t: 0 for t in BEGIN_TRANSITIONS}] + \
        [{t: 0 for t in MIDDLE_TRANSITIONS} for _ in range(match_states - 1)] + \
        [{t: 0 for t in END_TRANSITIONS}]

    for i, transitions in enumerate(transition_matrix):
        counts = transition_counts[i]
        totals = get_transition_totals(counts)
        possible_transitions = get_transition_possibilities(transitions)
        for key in transitions:
            transitions[key] = np.log((counts[key] + pseudocount) /
                                      (totals[key[0]] + (possible_transitions[key[0]] * pseudocount)))

    return transition_matrix


def get_transition_totals(transition_counts):
    """
    Given a dictionary of transition counts, returns a dictionary with the total number of times each type of transition from a particular state appears in the dictionary.

    Example:
        >>> get_transition_totals({'BM': 1, 'BI': 2, 'BD': 3})
        {'B': 6}
        >>> get_transition_totals({'BM': 1, 'BI': 2, 'BD': 3, 'MM': 4, 'MI': 5, 'MD': 6})
        {'B': 6, 'M': 15}
    """
    combined = {}

    for key, value in transition_counts.items():
        first_letter = key[0]

        if first_letter in combined:
            combined[first_letter] += value
        else:
            combined[first_letter] = value

    return combined


def get_transition_possibilities(transition_counts):
    """
    Given a dictionary of transition counts, returns a dictionary with the number of possible transitions from a particular state.

    Example:
        >>> get_transition_possibilities({'BM': 1, 'BI': 2, 'BD': 3})
        {'B': 3}
    """
    counts = {}

    for key in transition_counts.keys():
        first_letter = key[0]

        if first_letter in counts:
            counts[first_letter] += 1
        else:
            counts[first_letter] = 1

    return counts


def get_emission_counts(multiple_alignment):
    """
    Given a multiple alignment, returns a list of dictionaries, where each dictionary contains the number of times each amino acid appears in each column corresponding to a match state.
    """
    match_states, match_pattern = get_match_states(multiple_alignment)
    alignment_columns = get_alignment_columns(multiple_alignment)
    emission_counts_list = []

    for col, match in zip(alignment_columns, match_pattern):
        emission_counts = {a: 0 for a in AMINO_ACIDS}
        if match:
            for char in col:
                if char != '-':
                    emission_counts[char] += 1
            emission_counts_list.append(emission_counts)

    return emission_counts_list


def get_emission_matrix(multiple_alignment, pseudocount=0.1):
    """
    Given a multiple alignment, returns a list of dictionaries, where each dictionary contains the emission probabilities for each column corresponding to a match state.
    """
    emission_counts = get_emission_counts(multiple_alignment)
    emission_totals = [sum(emission_counts[i].values())
                       for i in range(len(emission_counts))]
    emission_matrix = []

    for j, col in enumerate(emission_counts):
        emission_probabilities = {
            a: np.log((col[a] + pseudocount) / (emission_totals[j] + len(AMINO_ACIDS * pseudocount))) for a in AMINO_ACIDS}
        emission_matrix.append(emission_probabilities)

    return emission_matrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment_file', type=str,
                        help='Path to the alignment file')
    args = parser.parse_args()
    sequences = ra.read_alignment(
        args.alignment_file, args.alignment_file[args.alignment_file.find('.') + 1:])

    # # Print transition matrix
    # print(f"\n{'=-' * 20}\n\nTRANSITION MATRIX\n\n{'=-' * 20}\n")
    # tm = get_transition_matrix(sequences, 1)
    # [print(f"{tm[i]}\n") for i in range(len(tm))]

    # # Print emission matrix
    # print(f"\n{'=-' * 20}\n\nEMISSION MATRIX\n\n{'=-' * 20}\n")
    # em = get_emission_matrix(sequences, 1)
    # [print(f"{em[i]}\n") for i in range(len(em))]

    hmm = ProfileHMM(sequences, 1)
