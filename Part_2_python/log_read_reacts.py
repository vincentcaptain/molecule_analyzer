'''
Owner: Qian Yang 06/20/2016
Last edit: Enze Chen 10/22/2018

This script extracts data from findmolecules.out and writes to the following files:
    reactperframe - rows are frames, columns are reaction IDs, data is number (SPARSE)
    reactdict - reactions are IDed in order of appearance
    reactbasis - rows are reactions, columns are molecule IDs, data is stoichiometry (SPARSE)
    reactantsbasis - rows are reactions, columns are molecules IDs, data is stoichiometry (SPARSE)
'''

import re
import operator

def getReacts(datafiles, maxframenum, moleculedict, newmoldict, bonded_reaction, bond_dist, filetags=[]):
    
    # Common dictionary over all files
    reactionbasisdict = {}  # Index for all reactions in order of appearance
    reactionbasis = {}  # Stoichiometry of all molecules in each reaction. Stored as sparse matrix
    reactantsbasis = {} # Stoichiometry of all reactants in each reaction. Stored as sparse matrix
    newreactionbasisdict = {}
    newreactionbasis = {}  # Stoichiometry of all molecules in each reaction. Stored as sparse matrix
    newreactantsbasis = {} # Stoichiometry of all reactants in each reaction. Stored as sparse matrix
    reactnum = 0
    num_atoms_involved = {}


    num_reaction = -1

    filecount = 0
    for datafile in datafiles:
        # ID purposes
        filecount += 1
        if len(filetags) < filecount:
            filetag = str(filecount)
        else:
            filetag = filetags(filecount - 1)



        reactlistperframe = {} # Index of reaction in each frame. Dictionary to form sparse list
        newreactlistperframe = {} # Index of reaction in each frame. Dictionary to form sparse list
        newdistancereactlistperframe = {}
        framenum = 1; # initialize framenum
        framenumhelper = -1; # initialize helper
        
        with open(datafile) as infile:
            for curr in infile:
                if re.match('Beginning frame', curr):
                    framenum +=1
                    # framenumhelper += 1
                    # if framenumhelper % 2:
                    #     framenum += 1
                        
                    if framenum > maxframenum:
                        print(framenum)
                        break

                    # 0-pad the matrix in case no reactions appear in last frame(s)
                    if framenum > maxframenum - 100:
                        key = str(framenum) + ',1'
                        reactlistperframe[key] = 0

                    if framenum % 6000 == 0:
                        print(curr.split(' '))
                elif re.match('Change bond: ', curr):
                    num_reaction +=1

                elif re.match('Local reaction', curr):
                    splitline = curr.split(': ')
                    splitline2 = splitline[1].split('\n')
                    reactionname = splitline2[0]

                elif re.match('Local bond reaction: ', curr):

                    if reactionname not in reactionbasisdict:
                        reactionbasisdict[reactionname] = len(reactionbasisdict)
                        reactnum += 1
                        
                        # separate left and right sides
                        [reactionLeft, reactionRight] = reactionname.split('=> ')
                        reactionLeftMols = reactionLeft.split('+ ')
                        reactionRightMols = reactionRight.split('+ ')
                
                        for mol in reactionLeftMols:
                            key = str(reactnum) + ',' + str(moleculedict[mol] + 1)
                            if(key in reactionbasis):
                                reactionbasis[key] -= 1
                                reactantsbasis[key] += 1
                            else:
                                reactionbasis[key] = -1
                                reactantsbasis[key] = 1
                                
                        for mol in reactionRightMols:
                            key = str(reactnum) + ',' + str(moleculedict[mol] + 1)
                            if(key in reactionbasis):
                                reactionbasis[key] += 1
                            else:
                                reactionbasis[key] = 1


                        splitline = curr.split(': ')
                        [Left_bond, Right_bond] = splitline[1].split('=> ')
                        Left_bond = Left_bond.split('+ ')
                        Right_bond = Right_bond.split('\n')
                        Right_bond = Right_bond[0].split('+ ')
                        new_Left = getBondedAtom(reactionLeftMols, Left_bond)
                        new_Right = getBondedAtom(reactionRightMols, Right_bond)
                        newreactionname = ''
                        for i in range(len(new_Left)):
                            mol = new_Left[i]
                            newreactionname += str(mol)
                            if i != len(new_Left)-1:
                                newreactionname += '+ '
                            key = str(reactnum) + ',' + str(newmoldict[mol] + 1)
                            if (key in newreactionbasis):
                                newreactionbasis[key] -= 1
                                newreactantsbasis[key] += 1
                            else:
                                newreactionbasis[key] = -1
                                newreactantsbasis[key] = 1

                        newreactionname += '=> '

                        for i in range(len(new_Right)):
                            mol = new_Right[i]
                            newreactionname += str(mol)
                            if i != len(new_Right)-1:
                                newreactionname += '+ '
                            key = str(reactnum) + ',' + str(newmoldict[mol] + 1)
                            if (key in newreactionbasis):
                                newreactionbasis[key] += 1
                            else:
                                newreactionbasis[key] = 1

                        num_atoms_involved[len(newreactionbasisdict)] = len(Right_bond)

                        newreactionbasisdict[newreactionname] = len(newreactionbasisdict)

                    try:
                        key = str(framenum) + ',' + str(reactionbasisdict[reactionname] + 1)
                        if (key in newreactlistperframe):
                            newreactlistperframe[key] += 1
                        else:
                            newreactlistperframe[key] = 1

                        # key = str(framenum) + ',' + str(reactionbasisdict[reactionname] + 1) + ',' + str(bonded_reaction[num_reaction])
                        # if (key in newdistancereactlistperframe):
                        #     newdistancereactlistperframe[key] += 1
                        # else:
                        #     newdistancereactlistperframe[key] = 1


                        key = str(framenum) + ',' + str(reactionbasisdict[reactionname] + 1)
                        if(key in reactlistperframe):
                            reactlistperframe[key] += 1
                        else:
                            reactlistperframe[key] = 1
                    except:
                        print('filecount: {}'.format(filecount))
                        print('framenum: {}'.format(framenum))
                        print('reactdictnum: '.format(reactionbasisdict[reactionname]))
                        raise

        reactperdist = {}
        bond_dist2 = []
        reaction = []
        for i in bond_dist.keys():
            t = i.split('\n')
            t = t[0].split(' ')
            t.append(bond_dist[i])
            bond_dist2.append(t)
        for i in newreactionbasisdict.keys():
            t = i.split(' => ')
            if t[0][0] != '(':
                t = t[0].split(' + ')
                reaction.append(t)
            else:
                continue
        for i in range(len(reaction)):
            for j in range(len(bond_dist2)):
                if bond_dist2[j][0] == reaction[i][0] and bond_dist2[j][1] == reaction[i][1]:
                    if bond_dist2[j][2] in reactperdist.keys():
                        reactperdist[bond_dist2[j][2]] += bond_dist2[j][3]
                    else:
                        reactperdist[bond_dist2[j][2]] = bond_dist2[j][3]




        print(reactperdist)
        numreacts = len(reactionbasisdict) 
        print('Current number of distinct reactions observed = {}, filecount = {}'.format(numreacts, filecount))
        #DumpReactlistperframe(reactlistperframe, 'reactperframe_' + filetag + '.dat', numreacts)
        DumpReactlistperframe(newreactlistperframe, 'newreactperframe_' + filetag + '.dat', numreacts)
        DumpDistanceReactionRate(bonded_reaction, reactperdist, 'distancereactionrate_' + filetag + '.dat', numreacts)

    print('Total number of distinct reactions observed: {}'.format(numreacts))
    #DumpReactionbasisdict(reactionbasisdict, 'reactdict_all.txt')
    DumpReactionbasisdict(newreactionbasisdict, 'newreactdict_all.txt')
    #DumpReactionbasis(reactionbasis, 'reactbasis_all.dat')
    #DumpReactantsbasis(reactantsbasis, 'reactantsbasis_all.dat')
    DumpReactionbasis(newreactionbasis, 'newreactbasis_all.dat')
    DumpReactantsbasis(newreactantsbasis, 'newreactantsbasis_all.dat')
    #DumpNumAtomInvolved(num_atoms_involved, 'num_atoms_involved_all.dat')

    return
    
def DumpReactlistperframe(reactlistperframe, dumpreactperframe, numreacts):
    f = open(dumpreactperframe, 'w')
    for i in reactlistperframe.keys():
        [numrow, numcol] = i.split(',')
        f.write("%s %s  %s\n" % (numrow, numcol, reactlistperframe[i]))
    f.close()
    print('Done with reactions per frame.')
    return

def DumpDistanceReactionRate(bonded_reaction, reactperdist, dumpreactperframe, numreacts):
    f = open(dumpreactperframe, 'w')
    for i in reactperdist.keys():
        if int(i) in bonded_reaction[0]:
            for j in range(len(bonded_reaction[1])):
                if int(i) == bonded_reaction[0][j]:
                    index = j
            f.write("%s %s %s %s\n" % (i, reactperdist[i], bonded_reaction[1][index] , bonded_reaction[1][j]/reactperdist[i]))
        else:
            f.write(
                "%s %s %s %s\n" % (i, reactperdist[i], 0, 0))
    f.close()
    print('Done with reactions per frame.')
    return

def DumpReactionbasisdict(reactionbasisdict, dumpdict):
    f = open(dumpdict, 'w')
    for key,value in sorted(reactionbasisdict.items(), key = operator.itemgetter(1)):
        f.write("%s: %s\n" % (value + 1,key))
    f.close()
    print('Done with reaction basis dictionary.')
    return
        
def DumpReactionbasis(reactionbasis, dumpbasis):    
    f = open(dumpbasis, 'w')
    for i in reactionbasis.keys():
        [numrow, numcol] = i.split(',')
        f.write("%s %s  %s\n" % (numrow, numcol, reactionbasis[i]))
    f.close()
    print('Done with reaction basis file.')
    return

def DumpReactantsbasis(reactantsbasis, dumpreactants):
    f = open(dumpreactants, 'w')
    for i in reactantsbasis.keys():
        [numrow, numcol] = i.split(',')
        f.write("%s %s  %s\n" % (numrow, numcol, reactantsbasis[i]))
    f.close()
    print('Done with reactants basis file. (end)')
    return

def DumpNumAtomInvolved(reactionbasisdict, dumpdict):
    f = open(dumpdict, 'w')
    for key in reactionbasisdict.keys():
        f.write("%s %s\n" % (key+1,reactionbasisdict[key]))
    f.close()
    print('Done with reaction basis dictionary.')
    return

def getBondedAtom(Left, Left_Bond):
    length = [len(Left), len(Left_Bond)]
    unique = set(length)
    if (len(unique)!= 1):
        print(length)
        print('Error: Not same length for Left, Right, Left_Bond or Right_Bond')
        return
    reactant = []
    store = {}
    for i in range(len(Left)):
        split = Left_Bond[i].split(' ')
        for j in range(len(split)-1):
            split[j] = int(split[j])
        if split[0] == 0:
            reactant.append(Left[i])
        else:
            store[i+1] = []
            for k in range(len(split) -1):
                store[i+1].append(split[k])
    visited = {}
    for k in store.keys():
        visited[k] = 0
    for k in store.keys():
        string = ''
        string = StringBondedAtom(string, k, store, visited, Left)
        if string !='':
            reactant.append(string)
    return reactant

def StringBondedAtom(string,m, store, visited, atom_list):
    if visited[m] == 1:
        return string
    string += '(' + str(atom_list[m-1])
    visited[m] = 1
    for i in range(len(store[m])):
        string = StringBondedAtom(string, store[m][i], store, visited, atom_list)
    string += ')'
    return string