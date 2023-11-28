import argparse
import read_alignment as ra

TRANSITIONS = ("MM", "MI", "MD", "IM", "II", "ID", "DM", "DI", "DD")
END = "END"


class ProfileHMM:
    def __init__(self, alphabet, states, transition, emission):
        self.alphabet = alphabet
        self.states = states
        self.transition = transition
        self.emission = emission


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
    Same as get_alignment_columns, but replaces dashes in non-match alignment columns with periods.
    """
    match_states, match_pattern = get_match_states(multiple_alignment)
    alignment_columns = get_alignment_columns(multiple_alignment)

    new_alignment_columns = []
    for col, match in zip(alignment_columns, match_pattern):
        if match:
            new_alignment_columns.append(col)
        else:
            new_alignment_columns.append(col.replace('-', '$'))

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


def get_transition_counts(multiple_alignment, pseudocount=0):
    """
    Given a multiple alignment, returns a list of dictionaries, where each dictionary contains the number of times each type of transition appears each column corresponding to a match state.
    """
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

    # Separate match groups into transitions
    transitions = []
    for group in match_groups:
        transitions.append([''.join(e) for e in list(zip(*group))])

    print(transitions)

    transition_counts = [
        get_transition_counts_helper(t, pseudocount) for t in transitions]

    return transition_counts


def get_transition_counts_helper(transitions, pseudocount):
    """
    Given a list of transitions, returns a dictionary with the number of times each transition appears in the list.

    Args:
        - transitions: a list containing strings whereby each pair of characters represents a transition

    Returns:
        - transition_counts: a dictionary with the number of times each type of transition appears in the list
    """
    transition_counts = {t: pseudocount for t in TRANSITIONS}
    for string in transitions:
        clean_str = string.replace("$", "")
        str_len = len(clean_str)
        if str_len < 2:
            return END
        elif str_len == 2:
            key = ('D' if clean_str[0] == '-' else 'M') + \
                ('D' if clean_str[1] == '-' else 'M')
            transition_counts[key] += 1
        else:
            pairs = [clean_str[i:i+2] for i in range(len(string) - 1)]

            start_pair = pairs.pop(0)
            start_key = ('D' if start_pair[0] == '-' else 'M') + 'I'
            transition_counts[start_key] += 1

            end_pair = pairs.pop(-1)
            end_key = 'I' + ('D' if end_pair[1] == '-' else 'M')
            transition_counts[end_key] += 1

            for _ in pairs:
                transition_counts['II'] += 1
    return transition_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment_file', type=str,
                        help='Path to the alignment file')
    args = parser.parse_args()
    sequences = ra.read_alignment(
        args.alignment_file, args.alignment_file[args.alignment_file.find('.') + 1:])
    # print(sequences)
    # print(get_match_states(sequences))
    # print(get_alignment_columns(sequences))
    # print(get_alignment_columns_with_match_states(sequences))
    print(get_transition_counts(sequences, 0))
