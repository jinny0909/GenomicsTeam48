"""Burrows-Wheeler Alignment - Protein Search
This is a simple implementation of the BMA for protein sequencing."""


def precalculation(ref):
    # create variables
    n = len(ref)
    SA = list()
    bwt = list()
    rotation_list = list()
    rotation_list_rev = list()

    count = dict()
    O = dict()
    O_rev = dict()
    proteinIDs = dict()
    pro_num = 1
    tots = dict()
    tots2 = dict()

    # alphabet
    alphabet = []
    for q in range(n):
        if ref[q] not in alphabet and ref[q] != '$':
            alphabet.append(ref[q])
        if ref[q] == '&':
            proteinIDs[q] = pro_num
            pro_num += 1

    for letter in alphabet:
        count[letter] = 0
        O[letter] = list()
        O_rev[letter] = list()
        tots[letter] = 0
        tots2[letter] = 0

    reverse = ref[0:n - 1][::-1] + '$'

    for i in range(n):
        rot = ref[i:] + ref[0:i]
        suffix = [rot, i]
        rotation_list.append(suffix)

        rev_rot = reverse[i:] + reverse[0:i]
        rev_suf = [rev_rot, i]
        rotation_list_rev.append(rev_suf)

        # update count
        if ref[i] != '$':
            for char in alphabet:
                if ref[i] < char:
                    count[char] = count[char] + 1

    # sort
    rotation_list.sort()
    rotation_list_rev.sort()

    for s in rotation_list:
        SA.append(s[1])
        bwt.append(s[0][-1:])

        if s[0][-1:] != '$':
            tots[s[0][-1:]] += 1
        for char in tots.keys():
            O[char].append(tots[char])

    for j in rotation_list_rev:
        if j[0][-1:] != '$':
            tots2[j[0][-1:]] += 1
        for char in tots2.keys():
            O_rev[char].append(tots2[char])

    return SA, count, O, O_rev, n, alphabet, proteinIDs


def inexact(W, z, SA, n, C, O, O_rev, alpha):
    D = calculateD(W, n, C, O_rev)
    SA_indices = inex_recur(W, len(W) - 1, z, 1, n - 1, D, C, O, alpha)
    return [SA[x] for x in SA_indices], D


# calculates the D array for a query
def calculateD(W, n, C, O_rev):
    k = 1
    l = n - 1
    z = 0
    D = list()
    for i in range(len(W)):
        k = C[W[i]] + O_rev[W[i]][k - 1] + 1
        l = C[W[i]] + O_rev[W[i]][l]
        if k > l:  # NOT been found
            k = 1
            l = n - 1
            z = z + 1
        D.append(z)

    return D


# # other implementation
# def calculateD(W, n, C, O):
#     z = 0
#     j = 0
#     D = list()
#     for i in range(len(W) - 1):
#         if W[j:i]:
#             z += 1
#             j = i + 1
#         D.append(z) #TODO: not substring
#     return D

def check_D(D, i):
    if i < 0:
        return 0
    else:
        return D[i]

def check_O(O, char, i):
    if i < 0:
        return 0
    else:
        return O[char][i]

def inex_recur(W, i, z, k, l, D, C, O, alphabet):
    tempset = set()

    # too many differences have been encountered
    if z < check_D(D, i):
        return set()
    if i < 0:  # entire query matched
        for m in range(k, l + 1):  ## TODO: check
            tempset.add(m)
        return tempset

    # indices for bwt match
    I = set()
    I = I.union(inex_recur(W, i - 1, z - 1, k, l, D, C, O,
                           alphabet))  # move down query string. Insertion
    for char in alphabet:
        # find the SA interval
        K = C[char] + check_O(O, char, k - 1) + 1
        L = C[char] + check_O(O, char, l)
        if K <= L:  # if the substring found
            I = I.union(inex_recur(W, i, z - 1, K, L, D, C, O, alphabet))  # Deletion
            if char == W[i]:  # correct alignment
                I = I.union(inex_recur(W, i - 1, z, K, L, D, C, O, alphabet))
            else:  # decrement z due to mismatch
                I = I.union(inex_recur(W, i - 1, z - 1, K, L, D, C, O, alphabet))
    return I


def bwa(read, array, c, O, O_rev, length, a, proteins, mismatch):
    print("Read: \"%s\"\nMax Difference Threshold: %d\n" % (read, mismatch))
    matches, D = inexact(read, mismatch, array, length, c, O, O_rev, a)

    # check if match overlaps with terminator
    threshold = len(read) * 0.8
    next = 0
    output = set()
    for match in matches:
        for p in proteins.keys():
            if match < p:
                next = p
                break

        # if terminator is in the middle
        if match + threshold > next:
            # remove
            matches.remove(match)
        else:
            output.add(proteins[next])

    print("%d match(es)\n" % (len(output)))
    print("Matched proteins: %s" % output)

    return output
