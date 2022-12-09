# referred to official bowtie repository & other implemntation
# sources: 
# sources also cited in the written report 

# base quality is not accounted in this implementation 
SEED = 28 

# construct Fm index object
class FmIndex: 
    def __init__(self, SA, bwt, count, O, O_rev, n): 
        self.SA = SA
        self.bwt = bwt
        self.count = count 
        self.O = O 
        self.O_rev = O_rev
        self.n = n 

class AlignmentRange: 
    def __init__(self, k, l, pattern): 
        self.k = k 
        self.numRanges =1 
        self.l = l 
        if pattern: 
            self.pattern = pattern
        else: 
            self.pattern = None 
 
class Range:
    def __init__(self, k=None, l=None, pattern=None): 

        if not k and not l: 
            self.numRanges = 0 
            self.alignmentsR = []
        else: 
            self.numRanges = 1 
            self.alignmentsR = [AlignmentRange(k, l, pattern)]
        
def precalculation_bowtie(ref):
    n = len(ref)
    SA = list()
    bwt = list()
    rotation_list = list()
    rotation_list_rev = list()

    count = dict()
    O = dict()
    O_rev = dict()
    tots = dict()
    tots2 = dict()

    alphabet = []
    for q in range(n):
        if ref[q] not in alphabet and ref[q] != '$':
            alphabet.append(ref[q])

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

    fmIndex = FmIndex(SA, bwt, count, O, O_rev, n)
    return  fmIndex, alphabet


def joinRanges(r1, r2):
    if not r1 and not r2: 
        return None
    elif not r1: 
        return r2 
    elif not r2: 
        return r1 

    totalRange = r1.numRanges +  r2.numRanges 

    ranges = Range()
    ranges.numRanges = totalRange
    ranges.alignmentsR = [0] * totalRange
    
    for i in range(r1.numRanges): 
        ranges.alignmentsR[i] = r1.alignmentsR[i]
    for i in range(r2.numRanges): 
        k = r2.alignmentsR[i].k
        l = r2.alignmentsR[i].l 
        for j in range(r1.numRanges): 
            if (k < ranges.alignmentsR[j].k): 
                if (l < ranges.alignmentsR[j].l): 
                    # delete jth item in alignmentsR for range
                    ranges.alignmentsR.pop(j)
                    ranges.alignmentsR[j] = r2.alignmentsR[i]
                    break 
                else: 
                    r2.alignmentsR.pop(i)
                    break
            elif (k < ranges.alignmentsR[j].l): 
                #should take quality into consideration, but will just ignore for the scope of this project 
                ranges.pop(r2.alignmentsR.pop(i))

            if ( j == r1.numRanges): 
                ranges.alignmentsR.append(r2.alignmentsR[i]) 
                ranges.numRanges += 1 

    return ranges 
                

def updateStartEnd(idx, pattern, length):
    i = len(pattern) - 1 
    c = pattern[i] 

    sp_aux = idx.count[c] - 1
    ep_aux = idx.count[c] - 1 
    i -= 1 

    while (sp_aux < ep_aux and i >= length): 
        c = pattern[i]
        sp_aux = idx.count[c] + idx.O[c][sp_aux] -  1
        ep_aux = idx.count[c] + idx.O[c][ep_aux] - 1 
        i -= 1 

    return sp_aux, ep_aux 

    
def inex_recur(W, i, minPos, z, k, l, C, O, alphabet, exact):

    if (i <= minPos): 
        if (z == 0 and exact == -1) or (z >= 0 and exact == 0) or (z < exact and exact > 0): 
            range = Range(k+1, l+1, W[len(W)-SEED:])
            return range 
        else: 
            return None 

    inter = Range()

    if (z > 0): 
        for char in alphabet:
            # find the SA interval
            K = C[char] + O[char][k] - 1
            L = C[char] + O[char][l] - 1 
            if K < L:  
                if char == W[i]:  # correct alignment
                    aux = inex_recur(W, i - 1, minPos, z, K, L, C, O, alphabet)
                else:  # z decrements because of mismatch 
                    aux = inex_recur(W, i - 1, minPos,  z - 1, K, L, C, O, alphabet, exact)

                inter = joinRanges(aux, inter)
    else:
        char = W[i]
        K = C[char] + O[char][k] - 1
        L = C[char] + O[char][l] -1
        if K < L: 
            inter = inex_recur(W, i-1, minPos,  z, K, L, C, O, alphabet, exact)


    return inter

reference = """KKKKLKKKKKKKKKKKKKKKKKKKKKKKKKKKKMFENITAAPADPILGLADLFRADERPGKINLGIGVYKDETGKTPVLTSVKKAEQYLLENETTKNYLGIDGIPEFGRCTQELLFGKGSALINDKRARTAQTPGGTGALRVAADFLAKNTSVKRVWVSNPSWPNHKSVFNSAGLEVREYAYYDAENHTLDFDALINSLNEAQAGDVVLFHGCCHNPTGIDPTLEQWQTLAQLSVEKGWLPLFDFAYQGFARGLEEDAEGLRAFAAMHKELIVASSYSKNFGLYNERVGACTLVAADSETVDRAFSQMKAAIRANYSNPPAHGASVVATILSNDALRAIWEQELTDMRQRIQRMRQLFVNTLQEKGANRDFSFIIKQNGMFSFSGLTKEQVLRLREEFGVYAVASGRVNVAGMTPDNMAPLCEAIVAVL"""
reference = reference + '$'
read = "MFENITAAPADPILGLADLFRADERPGKINLGIGVYKDETGKTPVLTSVKKAEQYLLENE"
rev_reference = reference[::-1]
rev_read = read[::-1]

# forward index 
forwardIdx, alphabet = precalculation_bowtie(reference)

# mirror index 
mirrorIdx, alphabet = precalculation_bowtie(rev_reference)

def exactMatch(fmi, pattern): 
    i = 0 
    sp, ep = updateStartEnd(fmi, pattern, 0)

    for i in range(sp+1, ep+1): 
        print("exact match result")
        print(fmi[i])


def case1(mirrorIdx, rev_read, mismatches, a):
    sp, ep = 0, int(len(rev_read) - SEED / 2)
    if (sp >= ep): 
        return 
    ranges = inex_recur(rev_read, int(len(rev_read) - SEED / 2 - 1), len(rev_read) - SEED - 1, mismatches, sp, ep, mirrorIdx.count , mirrorIdx.O, a, 0 )
    
    if not ranges: 
        return 

    pos = []
    for i in range(ranges.numRanges): 
        k = ranges.alignmentsR[i].k
        l = ranges.alignmentsR[i].l

        for j in range(k, l):
            
            posSuff = mirrorIdx.SA[j] - mirrorIdx.n + 1 
            pos.append(posSuff)

    return pos

def case2(forwardIdx, mirrorIdx, pattern, mismatch, a): 
    sp, ep = updateStartEnd(forwardIdx, pattern[0:SEED], int(SEED)/2)
    if(sp >= ep): 
        return 
    ranges = inex_recur(pattern[0:SEED], int(SEED/2 -1), -1, mismatch, sp, ep, forwardIdx.count, forwardIdx.O, a)
    if not ranges: 
        return 
    pos = []
    for i in range(ranges.numRanges): 
        seedPattern = ranges.alignmentsR[i].pattern 
        seedPattern = seedPattern[::-1]
        sp, ep = updateStartEnd(mirrorIdx, seedPattern, 0)
        for j in range(sp+1, ep+1):

            posForward = mirrorIdx.SA[j] - mirrorIdx.n +  1
            pos.append(posForward)

    return pos 

# when less than 1 mismatch allowed 

# set number of mismatch you want to allow 
mm = 1 
if mm == 0:
    case1(mirrorIdx, rev_read, 1, alphabet)
if mm == 1: 
    case1(mirrorIdx, rev_read, 1, alphabet)
    case2(forwardIdx, read, 1, alphabet)


