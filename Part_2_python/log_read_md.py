'''
Owner: Qian Yang 06/20/2016
Last edit: Enze Chen 10/22/2018

This script extracts data from molanal.out files and prints to the following files:
    molsperframe - rows are frames, data are total molecules in each frame
    molhistperframe - rows are frames, columns are molecule IDs (SPARSE)
    moleculedict - molecule full names are IDed in order of appearance
    moleculeatomdict - rows are molecules, columns are number of C,H in molecule
    bondsdict - rows are molecules, columns are number of CC, CH, HH bonds
'''

import re
import operator

def getMols(datafiles, maxframenum, filetags=[]):
    
    # Common dictionary over all files
    molmoleculedict = {} # Dictionary with all the molecules encountered to associate them an index
    newmoldict = {} # Dictionary with all the atom configurations encountered to associate them an index
    bond = {} # Dictionary which associate the index of the atom with all the indexes of the atoms it is bonded tot
    bond_prev = {} # Same as bond but for the previous frame
    configperframe = {} # Dictionary which associate the index of the atoms with their configurations
    config_prev = {} # Same as configperframe but for the previous frame
    newmolhistperframe = {} # For each time frame, count the number of atom of each configuration and store it
    newreactionbasisdict = {} # Dictionare with all the reactions between the atom configurations encountered to associate them an index
    reactnum = 0
    molmolhistperframe = {} # For each time frame, count the number of molecules and store it
    molmolhist = {} # Same as molmolhistperframe but only for the current frame
    newreactionbasis = {} # For each reaction, store by index which atom configurations participate in it and associate its stoechiometry
    newreactantsbasis = {} # Same as newreactionbasis but only for the reactants
    newreactlistperframe = {} # For each frame, store the number of time a reaction happened

    filecount = 0

    # Read the .out file line by line. This code goes twice through the file. The first allows to get some lists or
    # dictionaries which allow the second phase to be much faster.
    for datafile in datafiles:
        # ID purposes
        filecount += 1
        if len(filetags) < filecount:
            filetag = str(filecount)
        else:
            filetag = filetags(filecount - 1)

        framenum = 1;  # initialize framenum

        with open(datafile) as infile:
            for curr in infile:
                # At 'Beginning frame', dump molmolhistperframe to save memory and update framenum.
                if re.match('Beginning frame', curr):
                    if framenum == 1:
                        for i in molmolhist:
                            molmolhistperframe[str(framenum) + ',' + i] = molmolhist[i]
                        DumpMolhistperframe(molmolhistperframe, 'molmolhistperframe_' + filetag + '.dat', 0)
                        molmolhistperframe = {}

                    framenum += 1

                    if framenum > maxframenum:
                        # print(framenum)
                        break

                    if framenum % 6000 == 0:
                        print(curr.split(' '))
                    print(framenum)

                # At the beginning, store the molecules present
                elif re.match('Molecule: ', curr):
                    splitline = curr.split(': ')
                    splitline2 = splitline[1].split('\n')
                    molname = splitline2[0]
                    if molname not in molmoleculedict:
                        molmoleculedict[molname] = len(molmoleculedict)
                    key = str(framenum) + ',' + str(molmoleculedict[molname] + 1)
                    if (key in molmolhistperframe):
                        molmolhistperframe[key] += 1
                    else:
                        molmolhistperframe[key] = 1
                    key = str(molmoleculedict[molname] + 1)
                    if (key in molmolhist):
                        molmolhist[key] += 1
                    else:
                        molmolhist[key] = 1

                # Get the reactants and the products of the reaction with local features
                elif re.match('Local reaction', curr):
                    splitline = curr.split(': ')
                    splitline2 = splitline[1].split('\n')
                    [Left, Right] = splitline[1].split('=> ')
                    Left = Left.split('\n')
                    Right = Right.split('\n')
                    Left1 = Left[0].split('+ ')
                    Right1 = Right[0].split('+ ')

                # Using the reactants and products, and knowing if a bond is created or broken, the reaction is
                # stored in newreactbasisdict, if it is the first time it is encountered, and reactants and products are
                # stored in newreactionbasis and newreactantsbasis. Store the new atom configurations in newmoldict.
                elif re.match('Local bond reaction: ', curr):
                    splitline = curr.split(': ')
                    [Left_bond, Right_bond] = splitline[1].split('=> ')
                    Left_bond = Left_bond.split('+ ')
                    Right_bond = Right_bond.split('\n')
                    Right_bond = Right_bond[0].split('+ ')
                    new_Left = getBondedAtom(Left1, Left_bond)
                    new_Right = getBondedAtom(Right1, Right_bond)
                    for i in new_Left:
                        molname = i
                        if molname not in newmoldict:
                            newmoldict[molname] = len(newmoldict)
                    for i in new_Right:
                        molname = i
                        if molname not in newmoldict:
                            newmoldict[molname] = len(newmoldict)

                    newreactionname = ''
                    for i in range(len(new_Left)):
                        mol = new_Left[i]
                        newreactionname += str(mol)
                        if i != len(new_Left) - 1:
                            newreactionname += '+ '


                    newreactionname += '=> '

                    for i in range(len(new_Right)):
                        mol = new_Right[i]
                        newreactionname += str(mol)
                        if i != len(new_Right) - 1:
                            newreactionname += '+ '

                    if newreactionname not in newreactionbasisdict:
                        newreactionbasisdict[newreactionname] = len(newreactionbasisdict)
                        reactnum+=1

                        for i in range(len(new_Left)):
                            mol = new_Left[i]
                            key = str(reactnum) + ',' + str(newmoldict[mol] + 1)
                            if (key in newreactionbasis):
                                newreactionbasis[key] -= 1
                                newreactantsbasis[key] += 1
                            else:
                                newreactionbasis[key] = -1
                                newreactantsbasis[key] = 1
                        for i in range(len(new_Right)):
                            mol = new_Right[i]
                            key = str(reactnum) + ',' + str(newmoldict[mol] + 1)
                            if (key in newreactionbasis):
                                newreactionbasis[key] += 1
                            else:
                                newreactionbasis[key] = 1
                    key = str(framenum) + ',' + str(newreactionbasisdict[newreactionname] + 1)
                    if (key in newreactlistperframe):
                        newreactlistperframe[key] += 1
                    else:
                        newreactlistperframe[key] = 1

                # From the reaction with molecules, update molmolhist, by substracting reactants and adding products.
                elif re.match('Reaction with molecules', curr):
                    splitline = curr.split(': ')
                    [Left, Right] = splitline[1].split('=> ')
                    Left = Left.split('\n')
                    Right = Right.split('\n')
                    Left1 = Left[0].split('+ ')
                    Right1 = Right[0].split('+ ')
                    for i in Left1:
                        molname = i
                        key = str(molmoleculedict[molname] + 1)
                        molmolhist[key] -= 1
                    for i in Right1:
                        molname = i
                        if molname not in molmoleculedict:
                            molmoleculedict[molname] = len(molmoleculedict)
                        key = str(molmoleculedict[molname] + 1)
                        if (key in molmolhist):
                            molmolhist[key] += 1
                        else:
                            molmolhist[key] = 1

                elif re.match('End frame', curr):
                    for i in molmolhist:
                        molmolhistperframe[str(framenum) + ',' + i] = molmolhist[i]
                    DumpMolhistperframe(molmolhistperframe, 'molmolhistperframe_' + filetag + '.dat',1)
                    molmolhistperframe = {}

    filecount = 0

    num_atoms = 0

    #Going through a second time, but now it is easier with the info stored precedently.
    for datafile in datafiles:
        # ID purposes
        filecount += 1
        if len(filetags) < filecount:
            filetag = str(filecount)
        else:
            filetag = filetags(filecount - 1)

        framenum = 1; # initialize framenum


        with open(datafile) as infile:
            for curr in infile:
                if re.match('Beginning frame', curr):
                    if framenum ==1:
                        newmolhist = molPresent(bond_prev, config_prev, newmoldict)
                        for i in newmolhist:
                            newmolhistperframe[str(framenum) + ',' + str(i)] = int(newmolhist[i])
                    framenum += 1
                        
                    if framenum > maxframenum:
                        #print(framenum)
                        break
                    
                    if framenum % 6000 == 0:
                        print(curr.split(' '))
                    print(framenum)

                # At the first frame, store the configuration of each atom
                elif re.match('ID', curr):
                    splitline = curr.split(': ')
                    splitline2  = splitline[1].split('| ')
                    splitline3 = splitline2[1].split('\n')
                    molname = splitline3[0]
                    configperframe[str(framenum) + ',' + str(int(splitline2[0]) + 1)] = int(molname)
                    config_prev[int(splitline2[0]) + 1] = int(molname)
                    num_atoms += 1

                #At the first frame, store the atoms which are bonded to each atom
                elif re.match('Bond: ', curr):
                    splitline = curr.split(': ')
                    [Left, Right] = splitline[1].split(' is bonded to ')
                    Right = Right.split('\n')
                    Right = Right[0].split(' ')
                    Right.remove('')
                    string = []
                    for i in range(len(Right)):
                        Right[i] = int(Right[i])+1
                        string.append(Right[i])
                    bond[str(framenum) + ',' + str(int(Left) +1)] = string
                    bond_prev[str(int(Left)+1)] = Right

                #Modify the bonds of the atoms involved in the reaction
                elif re.match('Change bond: ', curr):
                    splitline = curr.split(': ')
                    [Left, Right] = splitline[1].split(' from ')
                    Left = Left.split(' and ')
                    Right = Right.split(' to ')
                    Right = Right[1].split(' \n')
                    Right = int(Right[0])
                    atom_involved = Left
                    if Right == 0:
                        bond_prev[str(int(Left[0])+1)].remove(int(Left[1])+1)
                        bond_prev[str(int(Left[1]) + 1)].remove(int(Left[0]) + 1)
                    else:
                        bond_prev[str(int(Left[0]) + 1)].append(int(Left[1]) + 1)
                        bond_prev[str(int(Left[1]) + 1)].append(int(Left[0]) + 1)

                # Modify the configuration of the atoms involved in the reaction
                elif re.match('Local reaction', curr):
                    splitline = curr.split(': ')
                    [Left, Right] = splitline[1].split('=> ')
                    Right = Right.split('\n')
                    Right1 = Right[0].split('+ ')
                    for i in range(len(Right1)):
                        config_prev[int(atom_involved[i]) +1] = int(Right1[i])

                # At the end of a frame, recalculate how many atom of each configuration of interest you have and store it.
                # Print only the first 2000 frames for bond and config.
                elif re.match('End frame', curr):
                    newmolhist = molPresent(bond_prev, config_prev, newmoldict)
                    for i in newmolhist:
                        newmolhistperframe[str(framenum) + ',' + str(i)] = int(newmolhist[i])
                    if framenum<2000:
                        for i in config_prev:
                            configperframe[str(framenum) + ',' + str(i)] = config_prev[i]
                        for i in bond_prev:
                            string = []
                            config = []
                            for j in range(len(bond_prev[i])):
                                config.append([config_prev[bond_prev[i][j]],bond_prev[i][j]])
                            config.sort()
                            for j in range(len(bond_prev[i])):
                                string.append(config[j][1])
                            bond[str(framenum) + ',' + i] = string
                    if framenum == 2000:
                        DumpBondperframe(bond, 'bondperframe_' + filetag + '.dat', 0)
                        DumpMolhistperframe(configperframe, 'configperframe_' + filetag + '.dat', 0)
                        bond = {}
                        configperframe = {}

    # Output all the necessary files
        DumpMolhistperframe(newmolhistperframe, 'newmolhistperframe_' + filetag + '.dat',0)
    DumpMoleculedict(newmoldict, 'newmoldict.txt')
    DumpMoleculedict(molmoleculedict, 'molmoleculedict_all.txt')


    numreacts = len(newreactionbasisdict)
    print('Current number of distinct reactions observed = {}, filecount = {}'.format(numreacts, filecount))
    DumpReactlistperframe(newreactlistperframe, 'newreactperframe_' + filetag + '.dat', numreacts)

    print('Total number of distinct reactions observed: {}'.format(numreacts))
    DumpReactionbasisdict(newreactionbasisdict, 'newreactdict_all.txt')
    DumpReactionbasis(newreactionbasis, 'newreactbasis_all.dat')
    DumpReactantsbasis(newreactantsbasis, 'newreactantsbasis_all.dat')

    return


def DumpMolhistperframe(molhistperframe, dumpmolhistfile,tag):
    if tag ==0:
        f = open(dumpmolhistfile, 'w')
        for i in molhistperframe.keys():
            if molhistperframe[i]!=0:
                [numrow, numcol] = i.split(',')
                f.write("%s %s %s\n" % (numrow, numcol, molhistperframe[i]))
        f.close()
        return
    else:
        f = open(dumpmolhistfile, 'a')
        for i in molhistperframe.keys():
            if molhistperframe[i] != 0:
                [numrow, numcol] = i.split(',')
                f.write("%s %s %s\n" % (numrow, numcol, molhistperframe[i]))
        f.close()
        return

def DumpBondperframe(molhistperframe, dumpmolhistfile, numspecies):
    f = open(dumpmolhistfile, 'w')
    for i in molhistperframe.keys():
        if molhistperframe[i]!=0:
            [numrow, numcol] = i.split(',')
            f.write("%s;%s;%s\n" % (numrow, numcol, molhistperframe[i]))
    f.close()
    print('Done with moleculehist per frame.')
    return

def DumpMoleculedict(moleculedict, dumpdictfile):
    f1 = open(dumpdictfile, 'w')
    for key,value in sorted(moleculedict.items(), key = operator.itemgetter(1)):
        f1.write("%s - %s  \n" % (value + 1, key))      # 1-index the molecule labels
    f1.close()
    print('Done with molecule and bonds dictionary.')
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

def are_bonded(atom_1, atom_2,father,bond):
    result = []
    for i in bond[str(atom_1)]:
        if i == atom_2:
            result.append(1)
        elif i in father:
            continue
        else:
            new_father = list(father)
            new_father.append(atom_1)
            a = are_bonded(i, atom_2, new_father,bond)
            if a !=0:
                result.append(a +1)
    if len(result)>0:
        return min(result)
    else:
        return 0

def unique(tab):
    result = []
    count = []
    for i in tab:
        if i in result:
            for j in range(len(result)):
                if result[j] == i:
                    count[j] +=1
        else:
            result.append(i)
            count.append(1)
    return [result, count]


def molPresent(bond, config, moldict):
    newmolhist = {}
    for i in moldict.keys():
        if i[0] != '(':
            for j in config.keys():
                if config[j] == int(i):
                    if moldict[i]+1 not in newmolhist:
                        newmolhist[moldict[i]+1] = 1
                    else:
                        newmolhist[moldict[i]+1] += 1
    for j in bond.keys():
        for k in bond[j]:
            if config[int(k)] > config[int(j)]:
                if '(' + str(config[int(j)]) + ' (' + str(config[int(k)]) +' ))' in moldict.keys():
                    if moldict['(' + str(config[int(j)]) + ' (' + str(config[int(k)]) +' ))']+1 in newmolhist.keys():
                        newmolhist[moldict['(' + str(config[int(j)]) + ' (' + str(config[int(k)]) +' ))']+1] += 0.5
                    else:
                        newmolhist[moldict['(' + str(config[int(j)]) + ' (' + str(config[int(k)]) +' ))']+1] = 0.5
            else:
                if '(' + str(config[int(k)]) + ' (' + str(config[int(j)]) +' ))' in moldict.keys():
                    if moldict['(' + str(config[int(k)]) + ' (' + str(config[int(j)]) +' ))']+1 in newmolhist.keys():
                        newmolhist[moldict['(' + str(config[int(k)]) + ' (' + str(config[int(j)]) +' ))']+1] += 0.5
                    else:
                        newmolhist[moldict['(' + str(config[int(k)]) + ' (' + str(config[int(j)]) +' ))']+1] = 0.5

    return newmolhist


def DumpReactlistperframe(reactlistperframe, dumpreactperframe, numreacts):
    f = open(dumpreactperframe, 'w')
    for i in reactlistperframe.keys():
        [numrow, numcol] = i.split(',')
        f.write("%s %s  %s\n" % (numrow, numcol, reactlistperframe[i]))
    f.close()
    print('Done with reactions per frame.')
    return

def DumpReactionbasisdict(reactionbasisdict, dumpdict):
    f = open(dumpdict, 'w')
    for key, value in sorted(reactionbasisdict.items(), key=operator.itemgetter(1)):
        f.write("%s: %s\n" % (value + 1, key))
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




