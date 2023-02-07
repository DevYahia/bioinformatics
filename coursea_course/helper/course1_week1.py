def pattern_count(text: str, pattern: str):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        sub_text = text[i: i + len(pattern)]
        if sub_text == pattern:
            count += 1
    return count


def frequency_table(text: str, k: int):
    freq_map = dict()
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i: i + k]
        if pattern in freq_map:
            freq_map[pattern] += 1
        else:
            freq_map[pattern] = 1
    return freq_map


def max_map(freq_map: dict):
    max_value = None
    for key, value in freq_map.items():
        if max_value is None:
            max_value = value
        elif value > max_value:
            max_value = value
    return max_value


def better_frequent_words(text: str, k: int):
    frequent_patterns = []
    freq_map = frequency_table(text, k)
    max_count = max_map(freq_map)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_count:
            frequent_patterns.append(pattern)
    return frequent_patterns


def reverse_complement(text: str):
    reversed_text = text[::-1]

    my_table = reversed_text.maketrans("ATGC", "TACG")

    c_dna = reversed_text.translate(my_table)

    return c_dna


def pattern_positions(genome: str, pattern: str):
    positions = []
    for i in range(len(genome) - len(pattern) + 1):
        sub_text = genome[i: i + len(pattern)]
        if sub_text == pattern:
            positions.append(i)
    return positions


def find_clumps(genome: str, k: int, L: int, t: int):
    patterns = []
    n = len(genome)
    for i in range(n - L + 1):
        window = genome[i: i + L]
        freq_map = frequency_table(window, k)
        for key, value in freq_map.items():
            if value >= t:
                patterns.append(key)
    return list(set(patterns))


def clumps_count(genome: str, k: int, L: int, t: int):
    count = 0
    kMers = []
    n = len(genome)
    for i in range(n - k + 1):
        kMers.append(genome[i: i + k])

    kMersIndices = dict()

    for i, kMer in enumerate(kMers):
        if kMer in kMersIndices:
            kMersIndices[kMer].append(i)
        else:
            kMersIndices[kMer] = [i]

    for kMer, indices in kMersIndices.items():
        indices_count = len(indices)
        if indices_count < t:
            continue
        for i in range(indices_count - t + 1):
            first_index = indices[i]
            last_index = indices[i + t - 1] + k
            diff = last_index - first_index
            if diff <= L:
                count += 1
                break

    return count

