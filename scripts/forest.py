#Written by Amanda Daigle
#Ernest Fraenkel's lab
#MIT Biological Engineering
#2013-2014


import sys, optparse, random, copy
import subprocess, tempfile
import networkx as nx
from operator import itemgetter



def score(value, mu, musquared):
    """
    Helper function for use in assigning negative prizes (when mu > 0)
    """
    if mu <= 0:
        raise ValueError('User variable mu is not greater than zero.')

    if value == 1:
        return 0
    else:
        newvalue = -float(value) #-math.log(float(value),2)
        if musquared: 
            newvalue = -(newvalue * newvalue)
        newvalue = newvalue * mu
        return newvalue
            
class PCSFInput(object):
    def __init__(self,prizeFile,edgeFile,confFile,dummyMode,knockout,garnet,shuffle,musquared,excludeT):
        """
        Converts input information into dictionaries to be used in the message passing algorithm
        
        INPUT: prizeFile - tab-delimited text file containing all proteins with prizes 
                           formatted like "ProteinName\tPrizeValue"
               edgeFile - tab-delimited text file containing edges in interactome and their weights
                          formatted like "ProteinA\tProteinB\tWeight\t
                          Directionality(U or D, optional)"
               confFile - text file containing values for all parameters. Should include the lines
                          "w=<value>", "D=<value>", and "b=<value>".
               dummyMode - a string that indicates which nodes in the interactome to connect the 
                           dummy node to. 'terminals'=connect to all terminals (default),  
                           'others'=connect to all nodes except for terminals, 'all'=connect to 
                           all nodes in the interactome, or a the path to a text file containing a 
                           list of proteins to connect to.
        
        OUTPUT: self.origPrizes - dictionary of proteins with prizes. {ProteinName: PrizeValue}
                self.negPrizes - dictionary of proteins with negative prizes calculated based on 
                                 degree. {ProteinName: PrizeValue}
                self.totalPrizes - dictionary of prizes as msgsteiner sees them (PrizeValue*beta)
                self.dirEdges - dictionary of dictionaries of all directed edges. 
                                {ProteinFROM: {ProteinTO: weight}}
                self.undirEdges - dictionary of dictionaries of all undirected edges. 
                                  {ProteinA: {ProteinB: weight}, ProteinB: {ProteinA: weight}}
                self.dummyNodeNeighbors - a list of all proteins that the dummy node should have
                                          edges to.
                self.interactomeNodes - a list of all nodes in the interactome
                self.w, self.b, self.D, self.n, self.mu, self.g, self.r, self.threads - parameters
        """
        if prizeFile==None or edgeFile==None:
            sys.exit('PCSF.py failed. Needs -p and -e arguments. Run PCSF.py -h for help.')
        warnings = 0
        selfedges = 0
        knockoutCount = 0
        #Check that dummyMode is a valid entry
        if dummyMode != 'terminals' and dummyMode != 'all' and dummyMode != 'others':
            try:
                dummyFile = open(dummyMode, 'rb')
            except:
                sys.exit('dummyMode value not recognized. Accepted values include "all", ' \
                         '"terminals", "others", or a path to a text file on your computer '\
                         'containing a list of proteins.')
        
        #Read configuration file to record parameters for this object
        print 'Reading text file containing parameters %s...' %confFile
        c = open(confFile, 'rb')
        for line in c:
            if line.startswith('w ='):
                w = line.strip().split()[-1]
            if line.startswith('b ='):
                b = line.strip().split()[-1]
            if line.startswith('D ='):
                D = line.strip().split()[-1]
            if line.startswith('mu ='):
                mu = line.strip().split()[-1]
            if line.startswith('n ='):
                n = line.strip().split()[-1]
            if line.startswith('r ='):
                r = line.strip().split()[-1]
            if line.startswith('g ='):
                g = line.strip().split()[-1]
            if line.startswith('threads = '):
                threads = line.strip().split()[1]
        c.close()
        try:
            mu = float(mu)
            #Warning for negative mu
            if mu < 0:
                print 'WARNING: mu = %f but must be a non-negative value. '\
                      'Changing mu to 0.\n' % mu
                warnings += 1
                mu = 0.0
        except:
            mu = 0.0
        try:
            n = float(n)
        except:
            n = 0.01 # Default n
        try:
            r = float(r)
        except:
            r = 0 # Default r
        try:
            g = float(g)
        except:
            g = 1e-3 # Default g
        try:
            threads = int(threads)
        except:
            threads = 1
        try:
            print 'Continuing with parameters w = %f, b = %f, D = %i, mu = %f, g = %f, n = %f.' \
                  %(float(w), float(b), int(D), mu, g, n)
        except:
            sys.exit('ERROR: There was a problem reading the file containing parameters. Please '\
                     'include appropriate values for w, b, D, and optionally mu, n, or g.')
                     
        self.w = float(w)
        self.b = float(b)
        self.D = int(D)
        self.mu = mu
        self.n = n
        self.r = r
        self.g = g
        self.threads = threads
        
        print 'Reading text file containing interactome edges: %s...' %edgeFile
        dirEdges = {}
        undirEdges = {}
        interactomeNodes = []
        try:
            e = open(edgeFile, 'rb')
        except IOError:
            sys.exit('ERROR: No such file %s, aborting program.\n' %edgeFile)
        line = e.readline()
        words = line.strip().split()
        #See if edgeFile contains directionality infomation in 4th column
        if len(words) == 3:
            col = 3
            print 'File does not contain direction information. Treating all edges as undirected\n'
        elif len(words) == 4:
            col = 4
            print 'File contains four columns. Fourth column will be interpreted as '\
                  'directionality information\n'
        else:
            print 'current line:', line
            sys.exit('ERROR: File containing interactome edges should have 3 or 4 columns: '\
                     'ProteinA\tProteinB\tWeight\tDirectionality(U or D). Protein names should '\
                     'not have spaces.')
        #to avoid printing out too many warnings, here is a way to tally up edited edges
        below_0,above1=0,0
        #Add each edge in edgeFile to undirEdges or dirEdges dictionary, appropriately
        while line:
            words = line.strip().split('\t')
            if len(words) != col:
                print 'current line:', line
                sys.exit('ERROR: All lines in the file containing the interactome edges should '\
                         'have the same number of columns. Protein names should not have spaces.')
            #Check for knockouts:
            if words[0] in knockout or words[1] in knockout:
                knockoutCount += 1
                line = e.readline()
                continue
            #Make sure edge weights are numbers between 0 and 0.99
            try:
                if float(words[2]) < 0:
                    below_0=below_0+1
                    words[2] = '0'
                    warnings += 1
                if float(words[2]) > 0.99:
                    above1=above1+1
                    words[2] = '0.99'
                    warnings += 1
            except:
                sys.exit('ERROR: Your interactome edge file include a non-numerical value in the '\
                         'third column. Aborting program.')
            #Check for self-edges
            if words[0] == words[1]:
                selfedges += 1
                line = e.readline()
                continue
            #Add all edges to undirEdges if there is no directionality information
            if col == 3:
                if words[0] not in undirEdges:
                    undirEdges[words[0]] = {}
                    undirEdges[words[0]][words[1]] = words[2]
                else:
                    if words[1] not in undirEdges[words[0]]:
                        undirEdges[words[0]][words[1]] = words[2]
                    else:
                        #If there are duplicate edges, keep the highest weight
                        if words[2] > undirEdges[words[0]][words[1]]:
                            undirEdges[words[0]][words[1]] = words[2]
                #Undirected edges represented as two directed edges
                if words[1] not in undirEdges:
                    undirEdges[words[1]] = {}
                    undirEdges[words[1]][words[0]] = words[2]
                else:
                    if words[0] not in undirEdges[words[1]]:
                        undirEdges[words[1]][words[0]] = words[2]
                    else:
                        #If there are duplicate edges, keep the highest weight
                        if words[2] > undirEdges[words[1]][words[0]]:
                            undirEdges[words[1]][words[0]] = words[2]
            #Add edges to dictionaries if there is directionality information
            if col == 4:
                if words[3] == 'U':
                    #If there are directed and undirected edges for the same protein pair, 
                    #only keep directed
                    if words[0] in dirEdges:
                        if words[1] in dirEdges[words[0]]:
                            line = e.readline()
                            continue
                    if words[0] not in undirEdges:
                        undirEdges[words[0]] = {}
                        undirEdges[words[0]][words[1]] = words[2]
                    else:
                        if words[1] not in undirEdges[words[0]]:
                            undirEdges[words[0]][words[1]] = words[2]
                        else:
                            #If there are duplicate edges, keep the highest weight
                            if words[2] > undirEdges[words[0]][words[1]]:
                                undirEdges[words[0]][words[1]] = words[2]
                    #Undirected edges represented as two directed edges
                    if words[1] not in undirEdges:
                        undirEdges[words[1]] = {}
                        undirEdges[words[1]][words[0]] = words[2]
                    else:
                        if words[0] not in undirEdges[words[1]]:
                            undirEdges[words[1]][words[0]] = words[2]
                        else:
                           #If there are duplicate edges, keep the highest weight
                            if words[2] > undirEdges[words[1]][words[0]]:
                                undirEdges[words[1]][words[0]] = words[2]
                elif words[3] == 'D':
                    if words[0] not in dirEdges:
                        dirEdges[words[0]] = {}
                        dirEdges[words[0]][words[1]] = words[2]
                    else:
                        if words[1] not in dirEdges[words[0]]:
                            dirEdges[words[0]][words[1]] = words[2]
                        else:
                            #If there are duplicate edges, keep the highest weight
                            if words[2] > dirEdges[words[0]][words[1]]:
                                dirEdges[words[0]][words[1]] = words[2]
                    #If there are directed and undirected edges for the same protein pair,
                    #keep only directed
                    if words[0] in undirEdges:
                        if words[1] in undirEdges[words[0]]:
                            del undirEdges[words[0]][words[1]]
                            del undirEdges[words[1]][words[0]]
                            if len(undirEdges[words[0]].keys()) == 0:
                                del undirEdges[words[0]]
                            if len(undirEdges[words[1]].keys()) == 0:
                                del undirEdges[words[1]]
                else:
                    print 'current line:', line
                    sys.exit('ERROR: The fourth column in the file containing interactome edges '\
                             'should only contain U or D.')
            #Keep track of nodes for dummyMode all and shufflePrizes
            if dummyMode =='all' or shuffle != 0:
                if words[0] not in interactomeNodes:
                    interactomeNodes.append(words[0])
                if words[1] not in interactomeNodes:
                    interactomeNodes.append(words[1])
            line = e.readline()
        e.close()
        self.interactomeNodes = interactomeNodes
        if above1 > 0:
            print 'WARNING!! All edgeweights should be a probability of protein '\
            'interaction. '+str(above1)+' of your edge weights include a number greater than 0.99.'\
            ' These were changed to 0.99...\n'
        if below_0 > 0:
            print 'WARNING!! All edgeweights should be a probability of protein '\
            'interaction. '+str(below_0)+' of your edge weights include a number below than 0. '\
            'These were changed to 0...\n'

        print 'Reading text file containing prizes: %s...\n' %prizeFile
        origPrizes = {}
        terminalTypes = {}
        try:
            p = open(prizeFile, 'rb')
        except IOError:
            sys.exit('ERROR: No such file %s, aborting program.\n' %prizeFile)
        #Count how many of these proteins are not in the interactome
        count = 0
        #Add each node in prizeFile to origPrizes dictionary
        line = p.readline()
        while line:
            words = line.strip().split()
            if len(words) != 2:
                print 'current line:', line
                sys.exit('ERROR: File containing prizes should have exactly two columns: '\
                         'ProteinName\tPrizeValue. Protein names should not have spaces.')
            #Increase count if this is not in the interactome
            if words[0] not in undirEdges and words[0] not in dirEdges:
                count += 1
            else:
                origPrizes[words[0]] = float(words[1])
                terminalTypes[words[0]] = 'Proteomic'
            line = p.readline()
        p.close()
        
        if garnet != None: 
            print 'Reading text file containing TF regression results: %s...\n' %garnet
            g = open(garnet, 'rb')
            line = g.readline()
            while line:
                words = line.strip().split()
                if len(words) != 2:
                    print 'current line:', line
                    sys.exit('ERROR: File containing TFs should have exactly two columns: '\
                         'TF_Name\tPrizeValue. TF names should not have spaces.')
                #Increase count if this is not in the interactome
                if words[0] not in undirEdges and words[0] not in dirEdges:
                    count += 1
                else:
                    prize = float(words[1])*self.n
                    #If the TF already has a prize value this will replace it.
                    origPrizes[words[0]] = prize
                    if words[0] in terminalTypes.keys():
                        terminalTypes[words[0]]+='_TF'
                    else:
                        terminalTypes[words[0]]='TF'
                line = g.readline()
            g.close()
        
        #Warning if supplied proteins were not in the interactome
        percentexcluded = (count/float(len(origPrizes.keys())+count)) * 100
        if percentexcluded > 90:
            sys.exit('ERROR: %i percent of your prize nodes are not included in the '\
                     'interactome! Make sure the protein names you are using are the same in your '\
                     'prize file as in your edge file. Aborting program.\n' %percentexcluded)
        elif percentexcluded > 0:
            print 'WARNING!! %.3f percent of your prize nodes are not included in the interactome!'\
                  ' These nodes were ignored. Make sure the protein names you are using are the '\
                  'same in your prize file as in your edge file. Continuing program...\n' \
                  %percentexcluded
            warnings += 1
        
        #Warning for self-edges
        if selfedges > 0:
            print 'WARNING: There were %i self-edges in your interactome. We ignored these '\
                  'edges.\n' %selfedges
            warnings += 1
        
        #Notice for knockouts
        if knockoutCount > 0:
            print 'There were %i edges connected to your knockout protein(s). We ignored these '\
                  'edges.\n' %knockoutCount
        
        print 'Input prize files and edge files have been successfully read.\n'
        
        #Connect dummy node to nodes in interactome, depending on dummyMode
        dummyNodeNeighbors = []
        if dummyMode == 'terminals':
            dummyNodeNeighbors = origPrizes.keys()
            print 'Dummy node has been added, with edges to all %i nodes which were assigned '\
                  'prizes.\n' %len(origPrizes.keys())
        elif dummyMode == 'all':
            dummyNodeNeighbors = interactomeNodes
            print 'Dummy node has been added, with edges to all %i nodes in the interactome.\n' \
                  %len(interactomeNodes)
        elif dummyMode == 'others':
            nonterminalNodes = []
            for node1 in undirEdges:
                for node2 in undirEdges[node1]:
                    if node2 not in origPrizes and node2 not in nonterminalNodes:
                        nonterminalNodes.append(node2)
                if node1 not in origPrizes and node1 not in nonterminalNodes:
                    nonterminalNodes.append(node1)
            for node1 in dirEdges:
                for node2 in dirEdges[node1]:
                    if node2 not in nonterminalNodes and node2 not in origPrizes:
                        nonterminalNodes.append(node2)
                if node1 not in nonterminalNodes and node1 not in origPrizes:
                    nonterminalNodes.append(node1)
            dummyNodeNeighbors = nonterminalNodes
            print 'Dummy node has been added, with edges to all %i nodes in the interactome '\
                  'which have not been assigned prizes.\n' %len(nonterminalNodes)
        else:
            #Keep track of how many genes on dummyNeighbors list are actually in interactome
            countNeighbors = 0.00
            numExcluded = 0
            line = dummyFile.readline()
            while line:
                line = line.strip()
                if line not in undirEdges and line not in dirEdges:
                    #protein not in interactome. Ignore edge but add to tally
                    numExcluded += 1
                else:
                    dummyNodeNeighbors.append(line)
                    countNeighbors += 1
                line = dummyFile.readline()
            dummyFile.close()
            if countNeighbors == 0:
                sys.exit('The file you provided for dummyMode does not contain any proteins in' \
                         'the interactome. Each line in your text file should contain the name of'\
                         'one protein. Make sure the names are the same as in the interactome.')
            percentexcluded = (numExcluded/countNeighbors) * 100
            #Warning if too many proteins in the file are excluded from dummyNodeNeighbors
            if percentexcluded > 0:
                print 'WARNING!! %i percent of the proteins listed in dummyNeighbors.txt are not '\
                      'included in the interactome! Make sure the protein names you are using '\
                      'are the same in this file as in the interactome. Continuing program...\n' \
                      %percentexcluded
                warnings += 1
            print 'Dummy node has been added, with edges to all %i nodes in the interactome '\
                  'listed in your dummyMode file.\n' %int(countNeighbors)

        self.terminalTypes = terminalTypes
        self.origPrizes = origPrizes
        self.dirEdges = dirEdges
        self.undirEdges = undirEdges
        self.dummyNodeNeighbors = dummyNodeNeighbors
        self.musquared = musquared
 
        self.assignNegPrizes(musquared, excludeT)
        
        if warnings > 0:
            print 'THERE WERE %s WARNING(S) WHEN READING THE INPUT FILES.\n' %warnings
        
    def assignNegPrizes(self, musquared, excludeT):     
        """        
        Scales original prizes by beta and adds negative prizes to penalize
        nodes with high degrees if mu > 0.
        Scales original prizes by beta if mu = 0.
        mu < 0 will cause score function to throw a ValueError.
        """
        negPrizes = {}
        totalPrizes = {}
        
        if self.mu != 0.0:
            print 'Adding negative prizes to nodes in interactome using mu parameter...'
            if musquared: print 'Negative prizes will be proportional to node degree^2.'
            if excludeT: print 'Terminals will retain their assigned prizes, no negative prizes.'
            DegreeDict = self.degreeNegPrize()
            for prot in self.origPrizes:
                if not excludeT:
                    try:
                        degree = DegreeDict[prot]
                        prize = (self.b * float(self.origPrizes[prot])) +\
                                 score(degree,self.mu,musquared)
                        negprize = score(degree,self.mu,musquared)
                        totalPrizes[prot] = prize
                        negPrizes[prot] = negprize
                    except KeyError:
                        continue
                else:
                    totalPrizes[prot] = self.b * float(self.origPrizes[prot])
                    negPrizes[prot] = 0
            for protein in DegreeDict:
                if protein not in self.origPrizes:
                    degree = DegreeDict[protein]
                    if degree > 0:
                        negprize = score(degree,self.mu,musquared)
                        if negprize != 0:
                            negPrizes[protein] = negprize
                            totalPrizes[protein] = negprize
        else:
            for prot in self.origPrizes:
                prize = float(self.origPrizes[prot]) * self.b
                totalPrizes[prot] = prize
                negPrizes[prot] = 0
       
        self.negPrizes = negPrizes
        self.totalPrizes = totalPrizes
    
    def degreeNegPrize(self):
        """
        Helper function for use in assigning negative prizes (when mu != 0)
        """
        G = nx.Graph()
        degreeDict = {}
        for protA in self.dirEdges:
            for protB in self.dirEdges[protA]:
                G.add_edge(protA, protB)
        for protA in self.undirEdges:
            for protB in self.undirEdges[protA]:
                G.add_edge(protA, protB)
        for node in G.nodes():
            degreeDict[node] = G.degree(node)
        return degreeDict

    def getInputInfo(self):
        """
        Prints the input information that this input object contains. Mostly used in debugging.
        """
        print 'The input prizes were', self.origPrizes
        print 'All undirected edges in the input interactome were', self.undirEdges
        print 'All directed edges in the input interactome were', self.dirEdges
        print 'The dummy node was connected to nodes: '+ str(self.dummyNodeNeighbors)
        print 'The parameters were: w= ' + self.w + ' b= ' +self.b+ ' D= ' + self.D + ' mu= '\
              + self.mu + ' r= ' + self.r + ' g= ' + self.g
    
    def runPCSF(self, msgpath, seed):
        """
        Passes the information in this input object to msgsteiner9, and returns the results.
        
        INPUT: msgpath - points to the directory where msgsteiner9 is held.
               
        RETURNS: edgeList - the contents of stdout from msgsteiner: a list of edges in the 
                            optimal Forest
                 info - the contents of stderr from msgsteiner: if all goes well, a report on
                             the optimization
        """
        print 'Preparing information to send to the message passing algorithm...\n'
        #Create a list of the input information for the msgsteiner subprocess
        input = tempfile.TemporaryFile()
        for edgeNode1 in self.dirEdges:
            for edgeNode2 in self.dirEdges[edgeNode1]:
                #directed edges are flipped so that they point towards the root node
                #Weights are converted to costs by using 1-weight 
                #(good for Psiquic, needs to be changed -log2(weight) for String)
                input.write('D %s %s %f\n' %(edgeNode2, edgeNode1, 
                            1-float(self.dirEdges[edgeNode1][edgeNode2])))
        edgesAdded = {}
        for edgeNode1 in self.undirEdges:
            for edgeNode2 in self.undirEdges[edgeNode1]:
                #keep track of undirected edges added so they are only included once in inputList
                try:
                    edgesAdded[edgeNode2][edgeNode1] = 2
                except KeyError:
                    input.write('E %s %s %f\n' %(edgeNode1, edgeNode2, 
                                1-float(self.undirEdges[edgeNode1][edgeNode2])))
                    edgesAdded[edgeNode1] = {edgeNode2: 1}
        for node in self.dummyNodeNeighbors:
            input.write('D %s DUMMY %.4f\n' %(node, self.w))
        for node in self.totalPrizes:
            input.write('W %s %f\n' %(node, float(self.totalPrizes[node])))
        input.write('W DUMMY 100.0\n')
        input.write('R DUMMY\n\n')

        input.seek(0)
        tempFileOut = open('tempFileData', 'w+')
        for line in input:
            tempFileOut.write(line)
        tempFileOut.close()
        print 'Input is processed. Piping to msgsteiner9 code...\n'
    
        #Run the msgsteiner code using Python's subprocess module
        try:
            with open(msgpath): pass
        except IOError:
            sys.exit('ERROR: The msgsteiner9 code was not found in the correct directory. '\
                     'Please use --msgpath to tell us the path to the msgsteiner9 code.' )
        

        #Run msgsteiner9 as subprocess. Using temporary files for stdin and stdout 
        #to avoid broken pipes when data is too big
        subprocArgs = [msgpath, '-d', str(self.D), '-t', '1000000', '-o', '-r', 
                       str(self.r), '-g', str(self.g), '-j', str(self.threads)]
        #Only supply seed to msgsteiner if one is given by user
        #Use the seed to set the -s (instance seed, which controls random noise on edge weights)
        #and the -z (message seed, which affects the message passing) msgsteiner seeds
        if seed != None:
            subprocArgs.append('-s')
            subprocArgs.append(str(seed))
            subprocArgs.append('-z')
            subprocArgs.append(str(seed))
        out = tempfile.TemporaryFile()
        input.seek(0) #return to first line of temporary file for reading
        subproc = subprocess.Popen(subprocArgs, bufsize=1, stdin=input, stdout=out, 
                                   stderr=subprocess.PIPE)
        errcode = subproc.wait()
        if errcode:
            errmess = subproc.stderr.read()
            sys.exit('ERROR: There was a problem running the message passing algorithm. <%s>: %s' \
                     %(errcode, errmess))
        print 'Message passing run finished with the parameters: w = %s, b = %s, D = %s, mu = %s' \
              ', g = %s\n' \
              %(self.w, self.b, self.D, self.mu, self.g)
        input.close()
        info = subproc.stderr.read()
        subproc.stderr.close()
        out.seek(0)
        edgeList = out.read()
        out.close()        
        return (edgeList, info)
        
class PCSFOutput(object):
    def __init__(self, inputObj, edgeList, info, outputpath, outputlabel, betweenness):
        """
        Takes the forest output given by msgsteiner and converts it to two networkx graphs.
        
        INPUT: inputObj - a reference to the PCSFInput object that stores the correct edges and 
                          prizes dictionaries for this output object
               edgeList - the edges in the forest output given by msgsteiner, contents of stdout
               info - stats about the msgsteiner run, contents of stderr
               outputpath - path to the directory where output files should be stored
               outputlabel - a label with which to name all of the output files for this run
               betweenness - a T/F flag indicating whether we should do the costly betweenness
                             calculation
               
        OUTPUT: self.optForest - a networkx digraph storing the forest returned by msgsteiner
                self.augForest - a networkx digraph storing the forest returned by msgsteiner, plus
                                 all of the edges in the interactome between nodes in that forest
                self.dumForest - a networkx digraph storing the dummy node edges in the optimal 
                                 forest
                self.inputObj - a reference to the PCSFInput object that created this output object
                <outputlabel>_info.txt - a text file containing the contents of stderr and info
        """
        #Write output stderr file before attempting to do anything else so the info is 
        #there if the program breaks
        err = open('%s/%s_info.txt' %(outputpath,outputlabel), 'wb')
        err.write(info)
        
        #Create networkx graph storing the result of msgsteiner
        optForest = nx.DiGraph()
        dumForest = nx.DiGraph()
        edges = edgeList.split('\n')
        for edge in edges:
            words = edge.split()
            if len(words) > 0:
                #If edge includes dummy node, only add to dumForest
                if words[0] == 'DUMMY':
                    sys.exit('ERROR: Tree returned from message passing algorithm has incorrect '\
                             'edges pointing towards the root dummy node!')
                if words[1] == 'DUMMY':
                    dumForest.add_edge(words[1], words[0])
                    continue
                #Add directed edges
                try:
                    optForest.add_edge(words[1], words[0], 
                                       weight=inputObj.dirEdges[words[1]][words[0]], 
                                       fracOptContaining=1.0)
                except KeyError:
                    #Add undirected edges
                    try:    
                        optForest.add_edge(words[1], words[0], 
                                           weight=inputObj.undirEdges[words[1]][words[0]], 
                                           fracOptContaining=1.0)
                        optForest.add_edge(words[0], words[1], 
                                           weight=inputObj.undirEdges[words[0]][words[1]], 
                                           fracOptContaining=1.0)
                    #edge not found in either dictionary
                    except KeyError:
                        sys.exit('ERROR: Edges were returned from the message passing algorithm '\
                                 'that were not found in the input data. Aborting program.')
        #Add prizes from input data as node attribute. Steiner nodes have prize zero.
        #Also add fracOptContaining as node attribute.
        terminalCount = 0
        for node in optForest.nodes():
            try:
                optForest.node[node]['prize'] = inputObj.totalPrizes[node]
                #Count terminal nodes (nodes that had prizes before mu)
                try:
                    orig = inputObj.origPrizes[node]
                    terminalCount += 1
                except KeyError:
                    pass
            except KeyError:
                optForest.node[node]['prize'] = 0
            optForest.node[node]['fracOptContaining'] = 1.0
            ##Added by sgosline: store terminal type to improve visualization
            try:
                ttype=inputObj.terminalTypes[node]
            except KeyError:
                ttype=''
            optForest.node[node]['TerminalType'] =ttype

        #Create networkx graph storing the "augmented forest", 
        #the result of msgsteiner plus all interactome edges between nodes present in the forest
        augForest = copy.deepcopy(optForest)
        for node in augForest.nodes():
            edges = {}
            try:
                edges.update(inputObj.undirEdges[node])
            except KeyError:
                pass
            try:
                edges.update(inputObj.dirEdges[node])
            except KeyError:
                pass
            for node2 in edges:
                if node2 in augForest.nodes():
                    if (node, node2) not in optForest.edges():
                        augForest.add_edge(node, node2, weight=edges[node2], fracOptContaining=0.0)
                    
        #Calculate betweenness centrality for all nodes in augmented forest
        if betweenness:
            betweenness = nx.betweenness_centrality(augForest)
            nx.set_node_attributes(augForest, 'betweenness', betweenness)
        else:
            for node in augForest.nodes():
                augForest.node[node]['betweenness'] = 0
        
        #Write info about results in info file
        err.write('\n')
        err.write('There were %i terminals in the interactome.\n' %len(inputObj.origPrizes.keys()))
        err.write('There are %i terminals in the optimal forest.\n' %terminalCount)
        err.write('The optimal forest has %i nodes total.\n' %len(optForest.nodes()))
        err.write('The roots in the optimal forest are ' + str([node for node in dumForest.nodes()
                  if node != 'DUMMY']) + '\n')
        err.write('Of these, ' + str([node for node in dumForest.nodes() if node != 'DUMMY' and 
                  node not in optForest.nodes()]) + ' are singletons and will not show up in '\
                  'optimalForest.sif since they have no leaves.\n')
        err.close()
        
        self.augForest = augForest
        self.optForest = optForest
        self.dumForest = dumForest
        self.inputObj = inputObj
            
    def writeCytoFiles(self, outputpath, outputlabel, cyto30):
        """
        Writes Cytoscape-supported files for viewing the output network.
        
        INPUT: outputpath - path to the directory where output files should be stored
               outputlabel - a label with which to name all of the output files for this run
        
        OUTPUT: <outputlabel>_optimalForest.sif - a file storing edges of optForest in Simple 
                                                  Interaction Format
                <outputlabel>_augmentedForest.sif - a file storing edges of augForest in Simple 
                                                    Interaction Format
                <outputlabel>_dummyForest.sif - a file storing edges of dumForest in Simple 
                                                Interaction Format
                <outputlabel>_nodeattributes.tsv - a tab-delimited file storing the node attributes
                <outputlabel>_edgeattributes.tsv - a tab-delimited file storing the edge attributes
        """
        
        if cyto30:
            #Write Simple Interaction Format files to store edges in a format supported by 
            #Cytoscape 3.0
            optSif= open('%s/%s_optimalForest.sif'%(outputpath,outputlabel), 'wb')
            augSif = open('%s/%s_augmentedForest.sif'%(outputpath,outputlabel), 'wb')
            dumSif = open('%s/%s_dummyForest.sif' %(outputpath, outputlabel), 'wb')
        
            #Write attribute files in a format supported by Cytoscape
            #The first lines of these files contains the variable names
            noa = open('%s/%s_nodeattributes.tsv'%(outputpath, outputlabel), 'wb')
            noa.write('Protein\tPrize\tBetweennessCentrality\t'\
                      'FractionOfOptimalForestsContaining\tTerminalType\n')
            eda = open('%s/%s_edgeattributes.tsv'%(outputpath, outputlabel), 'wb')
            eda.write('Edge\tWeight\tFractionOfOptimalForestsContaining\n')
        
            undirEdgesAdded = {}
            
            edgesSorted = self.augForest.edges(data=True)            
            edgesSorted.sort(key = itemgetter(0, 1))

            #iterate through edges to record edge types and edge attributes
            for (node1,node2,data) in edgesSorted:
                #Check if interaction between node1 and node2 is directed
                try:
                    w = self.inputObj.dirEdges[node1][node2]
                    augSif.write(node1+'\tpd\t'+node2+'\n')
                    eda.write(node1+' (pd) '+node2+'\t'+str(data['weight'])+'\t'+
                              str(data['fracOptContaining'])+'\n')
                    #Check if this edge is in optForest 
                    if data['fracOptContaining'] > 0.0:
                        optSif.write(node1+'\tpd\t'+node2+'\n')
                        dumSif.write(node1+'\tpd\t'+node2+'\n')
                #Else undirected
                except KeyError:
                    #Don't want to write undirected interactions twice (A pp B, B pp A).
                    try:
                        undirEdgesAdded[node2][node1] = 2
                    except KeyError:
                        augSif.write(node1+'\tpp\t'+node2+'\n')
                        eda.write(node1+' (pp) '+node2+'\t'+str(data['weight'])+'\t'+
                                  str(data['fracOptContaining'])+'\n')
                        if data['fracOptContaining'] > 0.0:
                            optSif.write(node1+'\tpp\t'+node2+'\n')
                            dumSif.write(node1+'\tpp\t'+node2+'\n')
                        if node1 in undirEdgesAdded:
                            undirEdgesAdded[node1][node2] = 1
                        else:
                            undirEdgesAdded[node1] = {node2:1}

            nodesSorted = self.augForest.nodes(data=True)
            nodesSorted.sort(key = itemgetter(0, 1))
            #iterate through nodes to record node attributes
            for (node,data) in nodesSorted:
                noa.write(node+'\t'+str(data['prize'])+'\t'+str(data['betweenness'])+'\t'+
                          str(data['fracOptContaining'])+'\t'+data['TerminalType']+'\n')

            dumSorted = self.dumForest.edges()
            dumSorted.sort(key = itemgetter(0, 1))        
            #Record dummy edges
            for (node1,node2) in dumSorted:
                if node1 == 'DUMMY':
                    dumSif.write(node1+'\tpd\t'+node2+'\n')
        
            optSif.close()
            augSif.close()
            dumSif.close()
            noa.close()
            eda.close()
            print 'Wrote output files for Cytoscape, in directory %s, with names starting with '\
                  '"%s".\n' %(outputpath, outputlabel)
                  
        else:
            #Write Simple Interaction Format files to store edges in a format supported by 
            #Cytoscape 2.8
            optSif= open('%s/%s_optimalForest.sif'%(outputpath,outputlabel), 'wb')
            augSif = open('%s/%s_augmentedForest.sif'%(outputpath,outputlabel), 'wb')
            dumSif = open('%s/%s_dummyForest.sif'%(outputpath,outputlabel), 'wb')
        
            #Write attribute files in a format supported by Cytoscape
            bcNoa = open('%s/%s_betweennessCentrality.noa'%(outputpath, outputlabel), 'wb')
            bcNoa.write('BetweennessCentrality (class=Double)\n')
            prizeNoa = open('%s/%s_prizes.noa'%(outputpath, outputlabel), 'wb')
            prizeNoa.write('Prize (class=Double)\n')
            fracNoa = open('%s/%s_fracOptContaining.noa'%(outputpath, outputlabel), 'wb')
            fracNoa.write('FractionOptimalForestsContaining (class=Double)\n')
            ttypeNoa =  open('%s/%s_termTypes.noa'%(outputpath, outputlabel), 'wb')
            ttypeNoa.write('TerminalType\n')
            
            #Write edge attribute files in a format supported by Cytoscape
            weightEda = open('%s/%s_weights.eda'%(outputpath, outputlabel), 'wb')
            weightEda.write('Weight (class=Double)\n')
            fracEda = open('%s/%s_fracOptContaining.eda'%(outputpath, outputlabel), 'wb')
            fracEda.write('FractionOptimalForestsContaining (class=Double)\n')
        
            undirEdgesAdded = {}
            edgesSorted = self.augForest.edges(data=True)            
            edgesSorted.sort(key = itemgetter(0, 1))            
            #iterate through edges to record edge types and edge attributes
            for (node1,node2,data) in edgesSorted:
                try:
                    #Check if interaction between node1 and node2 is directed
                    w = self.inputObj.dirEdges[node1][node2]
                    augSif.write(node1+'\tpd\t'+node2+'\n')
                    weightEda.write(node1+' (pd) '+node2+' = '+data['weight']+'\n')
                    fracEda.write(node1+' (pd) '+node2+' = '+str(data['fracOptContaining'])+'\n')
                    #Check if this edge is in optForest
                    if data['fracOptContaining'] == 1:
                        optSif.write(node1+'\tpd\t'+node2+'\n')
                        dumSif.write(node1+'\tpd\t'+node2+'\n')
                except KeyError:
                    #Don't want to write undirected interactions twice (A pp B, B pp A).
                    try:
                        undirEdgesAdded[node2][node1] = 2
                    except KeyError:
                        augSif.write(node1+'\tpp\t'+node2+'\n')
                        weightEda.write(node1+' (pp) '+node2+' = '+data['weight']+'\n')
                        fracEda.write(node1+' (pp) '+node2+' = '+str(data['fracOptContaining'])
                                      +'\n')
                        if data['fracOptContaining'] == 1:
                            optSif.write(node1+'\tpp\t'+node2+'\n')
                            dumSif.write(node1+'\tpp\t'+node2+'\n')
                        if node1 in undirEdgesAdded:
                            undirEdgesAdded[node1][node2] = 1
                        else:
                            undirEdgesAdded[node1] = {node2:1}
                    
            #iterate through nodes to record node attributes
            nodesSorted = self.augForest.nodes(data=True)
            nodesSorted.sort(key = itemgetter(0, 1))            
            for (node,data) in nodesSorted:
                bcNoa.write(node+' = '+str(data['betweenness'])+'\n')
                prizeNoa.write(node+' = '+str(data['prize'])+'\n')
                fracNoa.write(node+ ' = '+str(data['fracOptContaining'])+'\n')
                ttypeNoa.write(node+ ' = '+str(data['TerminalType'])+'\n')

            dumSorted = self.dumForest.edges()
            dumSorted.sort(key = itemgetter(0, 1))        
            #Record dummy edges
            for (node1,node2) in dumSorted:
                if node1 == 'DUMMY':
                    dumSif.write(node1+'\tpd\t'+node2+'\n')
        
            optSif.close()
            augSif.close()
            bcNoa.close()
            prizeNoa.close()
            fracNoa.close()
            weightEda.close()
            fracEda.close()
            print 'Wrote output files for Cytoscape, in directory %s, with names starting with '\
                  '"%s".\n' %(outputpath, outputlabel)
    
def mergeOutputs(PCSFOutputObj1, PCSFOutputObj2, betweenness, n1=1, n2=1):
    """
    Merges two PCSFOutput objects together. Creates a new PCSFOutput object whose graphs contain 
    all edges found in either original object, with updated fracOptContaining values and 
    betweenness values.
    
    INPUT: Two PCSFOutput objects, either individual objects or themselves results of merges.  
                Ideally, these output objects were created using the same interactome (though the 
                prizes or algorithm parameters may have been different). This is not enforced.
           betweenness - a T/F flag indicating whether to do the costly betweenness calculation
           n1,n2 - integers, the number of msgsteiner runs each PCSFOutput object is a result of 
                   (if one of PCSFOutputObj is the result of a merge, this should be >1).
                   
    RETURNS: A new PCSFOutput object, with all edges found in either original object, with updated 
                fracOptContaining and betweenness values. If a node or edge is found in both 
                original objects, the prize or weight in this object is copied from 
                PCSFOutputObj1. The inputObj reference in this object is the same as    
                PCSFOutputObj1.
    """
    print 'Merging outputs to give summary over %i algorithm runs...'%(n1+n2)
    mergedObj = copy.deepcopy(PCSFOutputObj1)
    #Update fracOptContaining for all edges in outputObj1
    for (node1,node2,data) in PCSFOutputObj1.optForest.edges(data=True):
        numRuns1 = data['fracOptContaining']*n1
        try:
            #if the edge is not in outputObj2 this will return a KeyError
            numRuns2 = PCSFOutputObj2.optForest[node1][node2]['fracOptContaining'] * n2
        except KeyError:
            numRuns2 = 0.0
        mergedObj.optForest[node1][node2]['fracOptContaining'] = (numRuns1 + numRuns2)/(n1+n2)
    #Update fracOptContaining for all nodes in outputObj1
    for (node, data) in PCSFOutputObj1.optForest.nodes(data=True):
        numRuns1 = data['fracOptContaining']*n1
        try:
            #if the node is not in outputObj2 this will return a KeyError
            numRuns2 = PCSFOutputObj2.optForest.node[node]['fracOptContaining'] * n2
        except KeyError:
            numRuns2 = 0.0
        mergedObj.optForest.node[node]['fracOptContaining'] = (numRuns1 + numRuns2)/(n1+n2)
    #Add optForest edges to mergedObj that appear in outputObj2 but not in outputObj1
    for (node1, node2, data) in PCSFOutputObj2.optForest.edges(data=True):
        try:
            dataM = mergedObj.optForest[node1][node2]
        except KeyError:
            numRuns2 = data['fracOptContaining'] * n2
            #If there are nodes in outputObj2 not included in 1, they will be added to 
            #mergedObj without error
            if node1 not in mergedObj.optForest.nodes():
                mergedObj.optForest.add_node(node1, 
                                             prize=PCSFOutputObj2.optForest.node[node1]['prize'],
                                             TerminalType=PCSFOutputObj2.optForest.node[node1]['TerminalType'],
                                             fracOptContaining=numRuns2/(n1+n2))
            if node2 not in mergedObj.optForest.nodes():
                mergedObj.optForest.add_node(node2, 
                                            prize=PCSFOutputObj2.optForest.node[node2]['prize'],
                                            TerminalType=PCSFOutputObj2.optForest.node[node2]['TerminalType'],
                                            fracOptContaining=numRuns2/(n1+n2))
            mergedObj.optForest.add_edge(node1, node2, weight=data['weight'], 
                                         fracOptContaining=numRuns2/(n1+n2))
    #Add dumForest edges to mergedObj that appear in outputObj2 but not in outputObj1
    for (node1, node2, data) in PCSFOutputObj2.dumForest.edges(data=True):
        try:
            dataM = mergedObj.dumForest[node1][node2]
        except KeyError:
            mergedObj.dumForest.add_edge(node1, node2)
    
    #Create augForest based on new optForest
    #Need to first copy optForest in case an edge previously included in augForest 
    #has a new fracOptContaining
    mergedObj.augForest = copy.deepcopy(mergedObj.optForest)
    for node in mergedObj.augForest.nodes():
        edges = {}
        try:
            edges.update(mergedObj.inputObj.undirEdges[node])
        except KeyError:
            pass
        try:
           edges.update(mergedObj.inputObj.dirEdges[node])
        except KeyError:
            #If a node found in mergedObj.optForest is not found in PCSFInputObj1's interactome, 
            #it is quietly ignored in making augForest
            pass
        for node2 in edges:
            if node2 in mergedObj.augForest.nodes():
                if (node, node2) not in mergedObj.optForest.edges():
                    mergedObj.augForest.add_edge(node, node2, weight=edges[node2], 
                                                 fracOptContaining=0.0)
                    
    #Calculate betweenness centrality for all nodes in augmented forest
    if betweenness:
        betweenness = nx.betweenness_centrality(mergedObj.augForest)
        nx.set_node_attributes(mergedObj.augForest, 'betweenness', betweenness)
        
    print 'Outputs were successfully merged.\n'
    return mergedObj
   
    
def shufflePrizes(PCSFInputObj, seed, excludeT):
    """
    Shuffles the prizes over all the nodes in PCSFInputObj.
    
    INPUT: a PCSFInput object
           seed - number to give to the random number generator
    RETURNS: a new PCSFInput object with the same prize values, but in randomly shuffled order. 
             All nodes will now have random prizes, but with the same distribution as the original 
             data.
    """
    print 'Prize values are being shuffled.\n'
    #Get number of original prizes
    numTruePrizes = len(PCSFInputObj.origPrizes)
    #Randomly choose which nodes will recieve prizes
    random.seed(seed)
    newNodes = random.sample(PCSFInputObj.totalPrizes.keys(), numTruePrizes)
    shuffledValues = dict(zip(newNodes,PCSFInputObj.origPrizes.values()))
    #Make a new PCSFInput object that contains all the same values as the original
    newPCSFInputObj = copy.deepcopy(PCSFInputObj)
    #Change the prizes to be the new dictionary
    newPCSFInputObj.origPrizes = shuffledValues
    newPCSFInputObj.assignNegPrizes(newPCSFInputObj.musquared,excludeT)
    return newPCSFInputObj
    
def noiseEdges(PCSFInputObj, seed, excludeT):
    """
    Adds gaussian noise to all edges in the PCSFInputObj prize dictionary.
    
    INPUT: a PCSFInput object
           seed - number to give to the random number generator
    RETURNS: a new PCSFInput object with with added gaussian noise to edge values
    """
    #Make a new PCSFInput object that contains all the same values as the original
    newPCSFInputObj = copy.deepcopy(PCSFInputObj)
    #Generate gaussian noise values, mean=0, stdev default=0.33 (edge values range between 0 and 1)
    random.seed(seed)
    if PCSFInputObj.n == None:
        dev = 0.333
    else:
        dev = PCSFInputObj.n
    for node1 in newPCSFInputObj.dirEdges:
        for node2 in newPCSFInputObj.dirEdges[node1]:
            newPCSFInputObj.dirEdges[node1][node2] = float(PCSFInputObj.dirEdges[node1][node2]) + \
                                                     random.gauss(0,dev)
    for node1 in newPCSFInputObj.undirEdges:
        for node2 in newPCSFInputObj.undirEdges[node1]:
            newPCSFInputObj.undirEdges[node1][node2] =float(PCSFInputObj.undirEdges[node1][node2])\
                                                      + random.gauss(0,dev)
    print 'Noise has been added to all edge values.\n'
    return newPCSFInputObj
    
def randomTerminals(PCSFInputObj, seed, excludeT):
    """
    Selects nodes with a similar degree distribution to the original terminals, and assigns the
    prizes to them.

    INPUT: a PCSFInput object
           seed - number to give to the random number generator
    RETURNS: a new PCSFInput object with a new list of terminals
    """
    #Only can do this if the interactome is big enough
    if len(PCSFInputObj.undirEdges) + len(PCSFInputObj.dirEdges) < 50:
        sys.exit("Cannot use --randomTerminals with such a small interactome.")
    #Make a new PCSFInput object that contains all the same values as the original but empty prizes
    newPCSFInputObj = copy.deepcopy(PCSFInputObj)
    newPCSFInputObj.origPrizes = {'':0}
    #degrees is a sorted list that will hold outdegree of every node in interactome
    degrees = []
    if len(PCSFInputObj.undirEdges) > 0 and len(PCSFInputObj.dirEdges) > 0:
        for node in PCSFInputObj.undirEdges:
            try:
                degrees.append((node,len(PCSFInputObj.undirEdges[node])+ \
                        len(PCSFInputObj.dirEdges[node])))
            except KeyError:
                degrees.append((node,len(PCSFInputObj.undirEdges[node])))
        for node in PCSFInputObj.dirEdges:
            if node not in PCSFInputObj.undirEdges:
                degrees.append((node,len(PCSFInputObj.dirEdges[node])))
    else:
        for node in PCSFInputObj.undirEdges:
            degrees.append((node,len(PCSFInputObj.undirEdges[node])))
        for node in PCSFInputObj.dirEdges:
            degrees.append((node,len(PCSFInputObj.dirEdges[node])))
    degrees.sort(key=itemgetter(1))
    #Find index of current terminal in degrees list
    for k,terminal in enumerate(PCSFInputObj.origPrizes):
        for i,value in enumerate(degrees):
            if terminal == value[0]:
                index = i
                break
        #Choose an index offset to select new terminal (distance from orig terminal in degrees list)
        #Make sure newly chosen terminal is not already chosen on a previous round
        newTerm = ''
        i = -1
        while newTerm in newPCSFInputObj.origPrizes and i<=10000:
            i+=1
            if seed != None:
                random.seed(seed+k+i)
            offset = int(random.gauss(0.0,100.0))
            newIndex = i + offset
            try:
                newNode = degrees[newIndex]
            except KeyError:
                #if offset points outside list, try loop again
                continue
            #To make truly random, need to choose randomly between all nodes with the same degree
            #Otherwise, ordering of dict iteration matters
            nodesWithSameDegree = []
            for node in degrees[newIndex:]:
                if node[1] == newNode[1]:
                    nodesWithSameDegree.append(node)
                else:
                    break
            for node in degrees[newIndex-1::-1]:
                if node[1] == newNode[1]:
                    nodesWithSameDegree.append(node)
                else:
                    break
            newTerm = random.choice(nodesWithSameDegree)[0]
        #if we've tried 10000 times, throw error to avoid infinite loop
        if newTerm in newPCSFInputObj.origPrizes:
            sys.exit('There was a problem with --randomTerminals. Aborting.')
        #Assign prize to newly chosen terminal
        newPCSFInputObj.origPrizes[newTerm] = PCSFInputObj.origPrizes[terminal]
    del newPCSFInputObj.origPrizes['']
    newPCSFInputObj.assignNegPrizes(newPCSFInputObj.musquared,excludeT)
    print 'New degree-matched terminals have been chosen.\n'
    return newPCSFInputObj

def changeValuesAndMergeResults(func, seed, inputObj, numRuns, msgpath, outputpath, outputlabel,
                                excludeT):
    """
    Changes the prizes/edges in the PCSFInput object according to func and runs the msgsteiner 
    algorithm, then merges the results together with the given PCSFOutput object. Writes 
    cytoscape files for final merged results.
    
    INPUT: func - the function which takes inputObj and changes the prize/edge values 
                  (i.e. shuffles or adds noise)
           inputObj - a PCSFInput object with original values, to be changed.
           numRums - the number of times to change the values and re-run msgsteiner
           msgpath - path to the directory where msgsteiner9 is kept
           outputpath - path to the directory where output files should be stored
           outputlabel - a label with which to name all of the output files for this run
           
    OUTPUT: <outputlabel>_changed_#_info.txt - a text file FOR EACH RUN containing the
                      contents of stderr for all msgsteiner runs
    RETURNS: merged - the PCSFOutput object that is a result of all the merges

    """ 
    print 'Preparing to change values %i times and get merged results of running the '\
          'algorithm on new values.\n' %numRuns
    i = 0
    while i <= numRuns-1:
        #Change prize/edge values and run msgsteiner with new values
        #NOTE: there will be an info file written for every run
        if seed != None:
            changedInputObj = func(inputObj, seed+i, excludeT)
            (newEdgeList, newInfo) = changedInputObj.runPCSF(msgpath, seed+i)
        else:
            changedInputObj = func(inputObj, seed, excludeT)
            (newEdgeList, newInfo) = changedInputObj.runPCSF(msgpath, seed)
        #By creating the output object with inputObj instead of changedInputObj, 
        #the prizes stored in the networkx graphs will be the ORIGINAL CORRECT prizes, 
        #not the changed prizes.
        if str(func)[10:23]  == 'shufflePrizes':
            changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath, 
                                          outputlabel+'_shuffledPrizes_%i'%i, 0)
        elif str(func)[10:20] == 'noiseEdges':
            changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath, 
                                          outputlabel+'_noisyEdges_%i'%i, 0)
        elif str(func)[10:25] == 'randomTerminals':
            changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath,
                                          outputlabel+'_randomTerminals_%i'%i, 0)
        if i == 0:
            #first run
            merged = changedOutputObj
        elif i == numRuns-1:
            #last run, merge results and calculate betweenness
            merged = mergeOutputs(merged, changedOutputObj, 1, i, 1)
        else:
            #Merge results of runs with the merged object containing all results so far
            merged = mergeOutputs(merged, changedOutputObj, 0, i, 1)
        i += 1
    #return merged outputobj
    return merged

def crossValidation(k, rep, PCSFInputObj, seed, msgpath, outputpath, outputlabel):
    """
    Seperates prizes into k "folds" and leaves those out of analysis. Reports what fraction of 
    held-out prize nodes were returned as steiner nodes.
    
    INPUT: k - the number of "folds" to seperate the prize nodes into.
           rep - Repetition of k-fold cv calculation we are currently on
           PCSFInputObj - a PCSF object with all prize nodes
           seed - number to give to the random number generator
           msgpath - path to msgsteiner code
           outputpath - path to the directory where output files should be stored
           outputlabel - a label with which to name all of the output files for this run
    
    OUTPUTS: File <outputlabel>_cvResults_<rep>.txt containing stats from the cv run 
             Files showing steiners and terminals for each of the intermediate solutions
    """ 
    print 'Running %i-fold cross validation (rep %i).\n'%(k,rep)
    prizes = PCSFInputObj.origPrizes.keys()
    if seed != None:
        random.seed(seed+rep)
    else: random.seed(None)
    random.shuffle(prizes)
    iterations = ''
    #Do k iterations
    for i in range(0,k):
        #File to store steiner and terminal nodes in results
        outputs = open('%s/%s_cvIntermediate_rep%ik%i.txt'%(outputpath,outputlabel,rep,i), 'wb')
        #select random prizes to hold out of this round
        hold_out = prizes[i:len(prizes):k]
        #keep track of these prizes which are returned in optimal network
        recovered = []
        #keep track of steiner nodes
        steiners = []
        #keep track of chosen terminals
        terminals = []
        newPCSFInputObj = copy.deepcopy(PCSFInputObj)
        for p in hold_out:
            #Remove held out original prize and update total prize to reflect only negPrize
            del newPCSFInputObj.origPrizes[p]
            newPCSFInputObj.totalPrizes[p] = newPCSFInputObj.negPrizes[p]
        (newEdgeList, newInfo) = newPCSFInputObj.runPCSF(msgpath, seed)
        #See if held out proteins appear in newEdgeList
        edges = newEdgeList.split('\n')
        for edge in edges:
            words = edge.split()
            if len(words) > 0:
                for node in (words[0], words[1]):
                    if node in hold_out:
                        if node not in recovered:
                            recovered.append(node)
                            steiners.append(node)
                    elif node not in prizes:
                            if node != 'DUMMY' and node not in steiners:
                                steiners.append(node)
                    else:
                        if node not in terminals:
                            terminals.append(node)
        #Write out lists for this fold's results
        outputs.write('Recovered Terminals\n')
        outputs.write(str(recovered))
        outputs.write('\nAll Steiner Nodes\n')
        outputs.write(str(steiners))
        outputs.write('\nTerminals\n')
        outputs.write(str(terminals))
        outputs.close
        #Return num of held-out terminals, num of recovered hold-outs, total num of Steiner nodes
        numRecovered = len(recovered)
        numSteiners = len(steiners)
        iterations = iterations + (str(i+1)+ '\t' + str(len(hold_out))+ '\t' + str(numRecovered)+\
                     '\t' + str(numSteiners) + '\n')
    results = open('%s/%s_cvResults_%i.txt'%(outputpath,outputlabel,rep), 'wb')
    results.write('Iteration\tNum of held-out terminals\tNum of recovered terminals\t'\
                  'Total num of Steiner nodes\n')
    results.write(iterations)
    results.close
        

def main():
    #Parsing arguments (run python PCSF.py -h to see all these decriptions)
    parser = optparse.OptionParser(description='Find multiple pathways within an interactome '\
        'that are altered in a particular condition using the Prize Collecting Steiner Forest '\
        'problem')
    #required arguments
    parser.add_option("-p", "--prize", dest='prizeFile', help='(Required) Path to the text file '\
        'containing the prizes. Should be a tab delimited file with lines: "ProteinName'\
        '\tPrizeValue"')
    parser.add_option("-e", "--edge", dest='edgeFile', help ='(Required) Path to the text file '\
        'containing the interactome edges. Should be a tab delimited file with 3 or 4 columns: '\
        '"ProteinA\tProteinB\tWeight(between 0 and 1)\tDirectionality(U or D, optional)"')
    #optional arguments
    parser.add_option("-c", "--conf", dest='confFile', help='Path to the text file containing '\
        'the parameters. Should be several lines that looks like: "ParameterName = '\
        'ParameterValue". Must contain values for w, b, D.  May contain values for optional '\
        'parameters mu, n, r, g. Default = "./conf.txt"', default='conf.txt')
    parser.add_option("-d","--dummyMode", dest='dummyMode', help='Tells the program which nodes '\
        'in the interactome to connect the dummy node to. "terminals"= connect to all terminals, '\
        '"others"= connect to all nodes except for terminals, "all"= connect to all '\
        'nodes in the interactome. If you wish you supply your own list of proteins, dummyMode '\
        'could also be the path to a text file containing a list of proteins (one per line). '\
        'Default = "terminals"', default='terminals')
    parser.add_option("--garnet", dest='garnet', help='Path to the text file containing '\
        'the output of the GARNET module regression. Should be a tab delimited file with 2 '\
        'columns: "TranscriptionFactorName\tScore". Default = "None"', default=None)
    parser.add_option("--musquared", action='store_true', dest='musquared', help='Flag to add '\
        'negative prizes to hub nodes proportional to their degree^2, rather than degree. Must '\
        'specify a positive mu in conf file.', default=False)
    parser.add_option("--excludeTerms", action='store_true', dest='excludeT', help='Flag to '\
        'exclude terminals when calculating negative prizes. Use if you want terminals to keep '\
        'exact assigned prize regardless of degree.', default=False)
    parser.add_option("--msgpath", dest='msgpath',  help='Full path to the message passing code. '\
        'Default = "<current directory>/msgsteiner9"', default='./msgsteiner9')
    parser.add_option("--outpath", dest = 'outputpath', help='Path to the directory which will '\
        'hold the output files. Default = this directory', default='.')
    parser.add_option("--outlabel", dest = 'outputlabel', help='A string to put at the beginning '\
        'of the names of files output by the program. Default = "result"', default='result')
    parser.add_option("--cyto30", action='store_true', dest='cyto30', help='Use this flag if '\
        'you want the output files to be amenable with Cytoscape v3.0 (this is the default).',
        default=True)
    parser.add_option("--cyto28", action='store_false', dest='cyto30', help='Use this flag if '\
        'you want the output files to be amenable with Cytoscape v2.8, rather than v3.0.')
    parser.add_option("--noisyEdges", dest='noiseNum', help='An integer specifying how many '\
        'times you would like to add noise to the given edge values and re-run the algorithm. '\
        'Results of these runs will be merged together and written in files with the word '\
        '"_noisyEdges_" added to their names. Default = 0', type='int', default=0)
    parser.add_option("--shuffledPrizes", dest='shuffleNum', help='An integer specifying how '\
        'many times you would like to shuffle around the given prizes and re-run the algorithm. '\
        'Results of these runs will be merged together and written in files with the word '\
        '"_shuffledPrizes_" added to their names. Default = 0', type='int', default=0)
    parser.add_option("--randomTerminals", dest='termNum', help='An integer specifying how many '\
        'times you would like to apply your given prizes to random nodes in the interactome (with'\
        ' a similar degree distribution) and re-run the algorithm. Results of these runs will be '\
        'merged together and written in files with the word "_randomTerminals_" added to their '\
        'names. Default = 0', type='int', default=0)
    parser.add_option("--knockout", dest='knockout', help='A list specifying protein(s) you '\
        'would like to "knock out" of the interactome to simulate a knockout experiment, '\
        "i.e. ['TP53'] or ['TP53', 'EGFR'].", type='string', default='[]')
    parser.add_option("-k", "--cv", dest='cv', help='An integer specifying the k value if you '\
        'would like to run k-fold cross validation on the prize proteins. Default = None.', \
        type='int', default=None)
    parser.add_option("--cv-reps", dest='cv_reps', help='An integer specifying how many runs of '\
        'cross-validation you would like to run. To use this option, you must also specify a -k '\
        'or --cv parameter. Default = None.', type='int', default=None)
    parser.add_option("-s", "--seed", dest='seed', help='An integer seed for the pseudo-random '\
        'number generators. If you want to reproduce exact results, supply the same seed. '\
        'Default = None.', type='int', default=None)
    
    (options, args) = parser.parse_args()
    
    #Check cv parameters do not conflict
    if options.cv_reps != None:
        if options.cv == None:
            sys.exit('You cannot use the --cv-reps option without also specifying a k parameter '\
                     'for k-fold cross validation.')

    #Process input, run msgsteiner, create output object, and write out results
    inputObj = PCSFInput(options.prizeFile,options.edgeFile, options.confFile, options.dummyMode,
                         options.knockout, options.garnet, options.shuffleNum,
                         options.musquared, options.excludeT)
    (edgeList, info) = inputObj.runPCSF(options.msgpath, options.seed)
    outputObj = PCSFOutput(inputObj,edgeList,info,options.outputpath,options.outputlabel,1)
    outputObj.writeCytoFiles(options.outputpath, options.outputlabel, options.cyto30)
    
    #Get merged results of adding noise to edge values
    if options.noiseNum > 0:
        merged = changeValuesAndMergeResults(noiseEdges, options.seed, inputObj, 
                                             options.noiseNum, options.msgpath, 
                                             options.outputpath, options.outputlabel,
                                             options.excludeT)
        merged.writeCytoFiles(options.outputpath, options.outputlabel+'_noisy', options.cyto30)
    
    #Get merged results of shuffling prizes
    if options.shuffleNum > 0:
        merged = changeValuesAndMergeResults(shufflePrizes, options.seed, inputObj, 
                                             options.shuffleNum, options.msgpath, 
                                             options.outputpath, options.outputlabel,
                                             options.excludeT)
        merged.writeCytoFiles(options.outputpath, options.outputlabel+'_shuffled', options.cyto30)

    #Get merged results of randomizing terminals
    if options.termNum > 0:
        merged = changeValuesAndMergeResults(randomTerminals,options.seed, inputObj,
                                             options.termNum, options.msgpath,
                                             options.outputpath, options.outputlabel,
                                             options.excludeT)
        merged.writeCytoFiles(options.outputpath, options.outputlabel+'_randomTerminals', 
                              options.cyto30)
    
    #If k is supplied, run k-fold cross validation
    if options.cv != None:
        if options.cv_reps == None:
            crossValidation(options.cv, 1, inputObj, options.seed, options.msgpath, \
                            options.outputpath, options.outputlabel)
        else:
            for i in range(0,options.cv_reps):
                crossValidation(options.cv, i+1, inputObj, options.seed, options.msgpath, \
                            options.outputpath, options.outputlabel)
    
if __name__ == '__main__':
    main()
