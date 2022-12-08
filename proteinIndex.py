import requests as r
from io import StringIO
from Bio import SeqIO  

def rotations(t):
    '''Return list of rotations of input string t'''
    tt = t * 2
    return [ tt[i:i+len(t)] for i in range(0, len(t)) ]

def bwm(t):
    '''Return lexicographically sorted list of t’s rotations'''
    return sorted(rotations(t))

def suffixArray(s):
    satups = sorted([(s[i:], i) for i in range(0, len(s))])
    arr = [] 
    for s in satups: 
      arr.append(s[1])
    del satups
    return arr

def bwtViaBwm(t):
    '''Given T, returns BWT(T) by way of the BWM'''
    return ''.join(map(lambda x: x[-1], bwm(t)))


class FMIndex:
  '''FM Index implementation using checkpointing and a suffix array sample'''

  #bwt_t: the BWT of the indexed string
  #cpStep: step for checkpoint indexes
  #checkpoints: array of dictionaries containing counts of each letter at each checkpointed index
  #F: representation of first column, maps letter -> range
  #saStep: step for SA offsets fo save
  #saSample: suffix array sample saving 1/saStep rows


  def __init__(self, t=None, cpStep=None, saStep=None, readFile=None):
    '''Build a new FM index'''
    if readFile != None:
      tots = self.loadFromFile(readFile)
      self.getF(tots)

    elif t != None and cpStep!= None and saStep != None and readFile == None:
      self.bwt_t = bwtViaBwm(t)
      self.cpStep = cpStep
      self.saStep = saStep
      tots = self.calcCheckpoints()
      self.getF(tots)
      self.calcSASample(t)
    else:
      print("ERROR: incorrect initialization")
    del tots

  
  def loadFromFile(self, proteinIndexFile):
    '''Load a previously saved FM index from file'''
    #Format: 
    #compressedBWT
    #!
    #cpStep
    #letter0:count0 letter1:count1....
    #(for all checkpoints)
    #!
    #saStep
    #idx self.saSample[idx]
    #(for call saSamples)
    tots = {}
    self.bwt_t = ""
    with open(proteinIndexFile, "r") as file:
      lines = file.readlines()

      bwtCompressed = lines[0]
      i = 0
      while(i < len(bwtCompressed)):
        if bwtCompressed[i] not in "1234567890": #character
          c = bwtCompressed[i]
          i+=1
          numStr = ""
          while(i < len(bwtCompressed) and bwtCompressed[i] in "1234567890"): #read digits
            numStr += bwtCompressed[i]
            i+=1
          if numStr != "" :
            self.bwt_t += c * int(numStr) #add run to bwt_t
            if c in tots.keys():
              tots[c] += int(numStr)
            else:
              tots[c] = int(numStr)

      self.cpStep = int(lines[2].strip()) #skip the !, begin reading checkpoints
      self.checkpoints = []
      j = 3
      while(lines[j][0]!="!"):
        checkpoint = {}
        if len(lines[j].strip()) > 0:
          lsplit = lines[j].strip().split(" ")
          for e in lsplit:
            char = e[0]
            num = int(str(e[1:]))
            checkpoint[char] = num
        j+=1
        self.checkpoints.append(checkpoint)

      j+=1 #skip the !
      self.saStep = int(lines[j].strip())
      j+=1
      self.saSample = {}
      for line in lines[j:]: #begin reading saSamples
        lsplit = line.strip().split(" ")
        self.saSample[int(lsplit[0])] = int(lsplit[1])
      file.close()
      return tots
  
  def calcSASample(self, t):
    '''get the suffix array and save every saStep samples'''
    satups = sorted([(t[i:], i) for i in range(0, len(t))])
    self.saSample = {} 
    rowi = 0
    for s in satups: 
      if s[1] % self.saStep == 0:
        self.saSample[rowi] = s[1]
      rowi += 1 
    del satups

  def calcCheckpoints(self):
    '''	calculate ranks for each letter and save checkpoints every cpStep'''
    tots	= {}
    self.checkpoints	=	[]
    for	i, c	in	enumerate(self.bwt_t):
      if i % self.cpStep == 0:
        self.checkpoints.append(tots.copy())
      if c in tots.keys():
        tots[c] += 1
      else:
        tots[c] = 1
    return tots
  
  def getF(self, tots):
    '''find first column ranges for each char'''
    self.F	=	{}
    totc	= 0
    for	c, count	in sorted(tots.items()):
      self.F[c]	=	(totc,	totc	+	count)
      totc	+=	count

  def resolveCountFromCheckpoint(self, c, idx):
    '''Resolve the count for character c at position idx in bwt_t using checkpoints'''
    cpStep = self.cpStep
    checkpoints = self.checkpoints
    bwt_t = self.bwt_t
    
    checki = round(idx/cpStep)
    if checki >= len(checkpoints):
      checki = len(checkpoints) - 1
    checkpointedIdx = checki * cpStep

    if c in checkpoints[checkpointedIdx].keys():
      count = checkpoints[checkpointedIdx][c]
    else:
      count = 0
    
    diff = 0
    for i in range(min(idx, checkpointedIdx), max(idx, checkpointedIdx)):
      if bwt_t[i] == c:
        diff += 1

    if idx < checkpointedIdx:
      count += diff
    else:
      count -= diff
    
    return count

  def LF(self, idx, c): 
    '''LF mapping using checkpoints'''
    return self.F[c][0] + self.resolveCountFromCheckpoint(c, idx) #get the first occurance of c in F. Get the rank of the occurance of c at idx by using checkpoints and add.

  def query(self, p):
    '''Find offsets of P in T using the FM index'''
    l = len(p)
    searchRange = (0, len(self.bwt_t))
    for i in range(l):
      c = p[l-1-i]
      lower = self.LF(searchRange[0], c)
      upper = self.LF(searchRange[1], c)
      searchRange = (lower, upper)
      if lower == upper: #not found since our range converged before we finished P
        return (-1, -1)
    
    return [self.resolveOffsets(i) for i in range(searchRange[0], searchRange[1])] #convert the range to offsets in original string

  def resolveOffsets(self, idx):
    '''Use the SA sample to resolve offset at position idx in bwt_t for the original string'''
    steps = 0
    i = idx 
    while self.bwt_t[i] != "$": #if we get all the way to $ we have reached the end of the string
      if i in self.saSample.keys():
        steps += self.saSample[i]
        break
      steps += 1
      i = self._lf(i, self.bwt_t[i])
    
    return steps
  
  def encode(self):
    '''Encode this FM index for easy write to file'''
    #Format: 
    #compressedBWT
    #!
    #cpStep
    #letter0count0 letter1count1....
    #(for all checkpoints)
    #!
    #saStep
    #idx self.saSample[idx]
    #(for call saSamples)

    lastC = self.bwt_t[0]
    count = 0
    bwtStr = ""
    for i, c in enumerate(self.bwt_t):
      if c == lastC:
        count+=1
      else:
        bwtStr += lastC + str(count)
        lastC = c
        count = 1

    bwtStr += lastC + str(count) + "\n"
    
    checkpointStr = str(self.cpStep) + "\n"
    for c in self.checkpoints:
      line = ""
      for letter in c.keys():
        line += str(letter) + str(c[letter]) + " "
      checkpointStr += line.strip() + "\n"
    
    saStr = str(self.saStep) + "\n"
    for idx in self.saSample.keys():
     saStr += str(idx) + " " + str(self.saSample[idx]) + "\n"
    
    return bwtStr +"!\n"+ checkpointStr +"!\n"+ saStr


def buildProteinIndex(proteinIDFile, cpStep, saStep):
  '''Build protein database by reading specified uniprot IDs from file and concatenating all sequences (using a terminator). Use to create an FM index'''
  pIDs = []
  with open(proteinIDFile, "r") as file:
    for line in file:
      pIDs.append(line.strip())
    file.close()

  baseUrl="http://www.uniprot.org/uniprot/"
  t = ""
  pos = 0
  proteinIDs = {}

  for pID in pIDs:
    currentUrl=baseUrl+pID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    Seq=StringIO(cData)
    pSeq=SeqIO.read(Seq,'fasta')
    proteinIDs[pos] = pSeq.id #map offset to the protein ID
    pos += len(pSeq.seq) + 1
    t += pSeq.seq + "&" #use & as terminator character for each protein sequence
  
  t += "$"
  del pIDs
  proteinIndex = FMIndex(t, cpStep, saStep)
  return proteinIndex

def saveProteinIndex(saveFile, proteinIndex):
  '''Write FM index to file'''
  with open(saveFile, "w") as file:
    file.write(proteinIndex.encode())
  file.close()
  print("Wrote to " + saveFile)
