import coursea_course.helper.constants as consts
import coursea_course.helper.course1_week1 as c1w1


def skew(genome: str):
    n = len(genome)
    total = [0]
    for i in range(n):
        if genome[i] == "G":
            total.append(total[i] + 1)
        elif genome[i] == "C":
            total.append(total[i] - 1)
        else:
            total.append(total[i])
    return total


def minimum_skew(genome: str):
    result = skew(genome)
    min_val = min(result)
    indices = []
    n = len(result)
    for i in range(n):
        if result[i] == min_val:
            indices.append(i)
    return indices


def hamming_distance(seq0: str, seq1: str):
    n = len(seq0)
    dist = 0
    for i in range(n):
        if seq0[i] != seq1[i]:
            dist += 1
    return dist


def approximate_pattern_matching(pattern: str, genome: str, d: int):
    n = len(genome)
    k = len(pattern)
    positions = []
    for i in range(n - k + 1):
        if hamming_distance(pattern, genome[i: i + k]) <= d:
            positions.append(i)
    return positions


def approximate_pattern_count(genome: str, pattern: str, d: int):
    matches = approximate_pattern_matching(pattern, genome, d)
    count = len(matches)
    return count


def immediate_neighbors(pattern: str):
    neighborhood = {pattern}
    n = len(pattern)
    for i in range(n):
        nucleotide = pattern[i]
        _nucleotides = consts.nucleotides.difference(nucleotide)
        for _n in _nucleotides:
            neighbor = pattern[:i] + _n + pattern[i + 1:]
            neighborhood.add(neighbor)
    return neighborhood


def neighbors(pattern: str, d: int):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return consts.nucleotides
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for neighbor in suffix_neighbors:
        if hamming_distance(pattern[1:], neighbor) < d:
            for nucleotide in consts.nucleotides:
                neighborhood.add(nucleotide + neighbor)
        else:
            neighborhood.add(pattern[0] + neighbor)
    return neighborhood


def neighbors_exact(pattern: str, d: int):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return consts.nucleotides
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for neighbor in suffix_neighbors:
        if hamming_distance(pattern[1:], neighbor) == d - 1:
            for nucleotide in consts.nucleotides:
                neighborhood.add(nucleotide + neighbor)
        elif hamming_distance(pattern[1:], neighbor) == d:
            neighborhood.add(pattern[0] + neighbor)
    return neighborhood


def frequent_words_with_mismatches(genome: str, k: int, d: int):
    patterns = []
    freq_map = dict()
    n = len(genome)
    for i in range(n - k + 1):
        pattern = genome[i: i + k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map[neighbor] += 1
            else:
                freq_map[neighbor] = 1
    m = c1w1.max_map(freq_map)
    for key, value in freq_map.items():
        if value == m:
            patterns.append(key)
    return patterns


def frequent_words_with_mismatches_and_reverse_complements(genome: str, k: int, d: int):
    patterns = []
    freq_map = dict()
    n = len(genome)
    for i in range(n - k + 1):
        pattern = genome[i: i + k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map[neighbor] += 1
            else:
                freq_map[neighbor] = 1
        pattern_rc = c1w1.reverse_complement(pattern)
        neighborhood = neighbors(pattern_rc, d)
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map[neighbor] += 1
            else:
                freq_map[neighbor] = 1
    m = c1w1.max_map(freq_map)
    for key, value in freq_map.items():
        if value == m:
            patterns.append(key)
    return patterns
