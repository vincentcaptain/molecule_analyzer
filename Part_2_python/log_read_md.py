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
    moleculedict = {}   # Index for all molecules in order of appearance, stored in dictionary
    # moleculeatomdict = []   # keys are moleculedict value, values are [numC, numH], arranged in array not dict
    # bondsdict = {}      # Index for the number of CC, CH, HH bonds (columns) per molecule (row)
    molmoleculedict = {}
    newmoldict = {}
    bond = {}
    bond_prev = {}
    configperframe = {}
    config_prev = {}
    newmolhistperframe = {}
    count_of_bond = {}
    bonded_reaction = []
    bond_dist = {}
    newreactionbasisdict = {}
    reactnum = 0
    molmolhistperframe = {}
    molmolhist = {}
    newreactionbasis = {}
    newreactantsbasis = {}
    newreactlistperframe = {}
    reaction_creat = []


    filecount = 0

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

                elif re.match('Local reaction', curr):
                    splitline = curr.split(': ')
                    splitline2 = splitline[1].split('\n')
                    reactionname = splitline2[0]
                    [Left, Right] = splitline[1].split('=> ')
                    Left = Left.split('\n')
                    Right = Right.split('\n')
                    Left1 = Left[0].split('+ ')
                    Right1 = Right[0].split('+ ')


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
                # try:
                    key = str(framenum) + ',' + str(newreactionbasisdict[newreactionname] + 1)
                    if (key in newreactlistperframe):
                        newreactlistperframe[key] += 1
                    else:
                        newreactlistperframe[key] = 1

                    # key = str(framenum) + ',' + str(reactionbasisdict[reactionname] + 1) + ',' + str(bonded_reaction[num_reaction])
                    # if (key in newdistancereactlistperframe):
                    #     newdistancereactlistperframe[key] += 1
                    # else:
                    #     newdistancereactlistperframe[key] = 1


                    # key = str(framenum) + ',' + str(reactionbasisdict[reactionname] + 1)
                    # if (key in reactlistperframe):
                    #     reactlistperframe[key] += 1
                    # else:
                    #     reactlistperframe[key] = 1
                # except:
                #     print('filecount: {}'.format(filecount))
                #     print('framenum: {}'.format(framenum))
                #     print('reactdictnum: '.format(reactionbasisdict[reactionname]))
                #     raise


                    # if reactionname not in reactionbasisdict:
                    #     reactionbasisdict[reactionname] = len(reactionbasisdict)
                    #     reactnum += 1


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

    for i in newreactionbasisdict.keys():
        t = i.split(' => ')
        if t[0][0] != '(':
            t = t[0].split(' + ')
            reaction_creat.append(t)
        else:
            continue

    filecount = 0

    num_atoms = 0
    for datafile in datafiles:
        # ID purposes
        filecount += 1
        if len(filetags) < filecount:
            filetag = str(filecount)
        else:
            filetag = filetags(filecount - 1)
        
        # Read data and parse each line
        # molsperframe = []   # Total number of molecules in each frame
        molhistperframe = {}    # Number of each molecule in each frame. Dictionary to form sparse list
        molhist = {}
        framenum = 1; # initialize framenum


        with open(datafile) as infile:
            for curr in infile:

                if re.match('Beginning frame', curr):
                    if framenum ==1:
                        # DumpMolhistperframe(molmolhistperframe, 'molmolhistperframe_' + filetag + '.dat',0)
                        # molmolhistperframe = {}
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


                elif re.match('ID', curr):
                    splitline = curr.split(': ')
                    splitline2  = splitline[1].split('| ')
                    splitline3 = splitline2[1].split('\n')
                    molname = splitline3[0]
                    # if molname not in moleculedict:
                    #     moleculedict[molname] = len(moleculedict)
                    # key = str(framenum) + ',' + str(moleculedict[molname] + 1)
                    # if(key in molhistperframe):
                    #     molhistperframe[key] += 1
                    # else:
                    #     molhistperframe[key] = 1
                    # key = str(moleculedict[molname] + 1)
                    # if (key in molhist):
                    #     molhist[key] += 1
                    # else:
                    #     molhist[key] = 1
                    configperframe[str(framenum) + ',' + str(int(splitline2[0]) + 1)] = int(molname)
                    config_prev[int(splitline2[0]) + 1] = int(molname)
                    num_atoms += 1
                    #REMOVE MOLHIST


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
                    #CALCULATE NEWMOLHIST


                # elif re.match('Molecule: ', curr):
                #     splitline = curr.split(': ')
                #     splitline2 = splitline[1].split('\n')
                #     molname = splitline2[0]
                #     if molname not in molmoleculedict:
                #         molmoleculedict[molname] = len(molmoleculedict)
                #     key = str(framenum) + ',' + str(molmoleculedict[molname] + 1)
                #     if (key in molmolhistperframe):
                #         molmolhistperframe[key] += 1
                #     else:
                #         molmolhistperframe[key] = 1
                #     key = str(molmoleculedict[molname] + 1)
                #     if (key in molmolhist):
                #         molmolhist[key] += 1
                #     else:
                #         molmolhist[key] = 1


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
                        #bonded_reaction.append(-1)
                    else:
                        #bonded_reaction.append(are_bonded(int(Left[0])+1, int(Left[1]) +1 , [], bond_prev))
                        bond_prev[str(int(Left[0]) + 1)].append(int(Left[1]) + 1)
                        bond_prev[str(int(Left[1]) + 1)].append(int(Left[0]) + 1)


                # elif re.match('Atoms involved: ', curr):
                #     splitline = curr.split(': ')
                #     splitline2 = splitline[1].split(' ')
                #     # atom_involved = splitline2

                elif re.match('Local reaction', curr):
                    splitline = curr.split(': ')
                    [Left, Right] = splitline[1].split('=> ')
                    # Left = Left.split('\n')
                    Right = Right.split('\n')
                    # Left1 = Left[0].split('+ ')
                    Right1 = Right[0].split('+ ')
                    # for i in Left1:
                    #     molname = i
                    #     key = str(moleculedict[molname] + 1)
                    #     molhist[key] -= 1
                    # for i in Right1:
                    #     molname = i
                    #     if molname not in moleculedict:
                    #         moleculedict[molname] = len(moleculedict)
                    #     key = str(moleculedict[molname] + 1)
                    #     if (key in molhist):
                    #         molhist[key] += 1
                    #     else:
                    #         molhist[key] = 1
                    for i in range(len(Right1)):
                        config_prev[int(atom_involved[i]) +1] = int(Right1[i])
                    #REMOVE PART ON MOLECULEDICT

                # elif re.match('Local bond reaction: ', curr):
                #     splitline = curr.split(': ')
                #     [Left_bond, Right_bond] = splitline[1].split('=> ')
                #     Left_bond = Left_bond.split('+ ')
                #     Right_bond = Right_bond.split('\n')
                #     Right_bond = Right_bond[0].split('+ ')
                #     new_Left = getBondedAtom(Left1, Left_bond)
                #     new_Right = getBondedAtom(Right1, Right_bond)
                #     for i in new_Left:
                #         molname = i
                #         if molname not in newmoldict:
                #             newmoldict[molname] = len(newmoldict)
                #     for i in new_Right:
                #         molname = i
                #         if molname not in newmoldict:
                #             newmoldict[molname] = len(newmoldict)


                # elif re.match('Reaction with molecules', curr):
                #     splitline = curr.split(': ')
                #     [Left, Right] = splitline[1].split('=> ')
                #     Left = Left.split('\n')
                #     Right = Right.split('\n')
                #     Left1 = Left[0].split('+ ')
                #     Right1 = Right[0].split('+ ')
                #     for i in Left1:
                #         molname = i
                #         key = str(molmoleculedict[molname] + 1)
                #         molmolhist[key] -= 1
                #     for i in Right1:
                #         molname = i
                #         if molname not in molmoleculedict:
                #             molmoleculedict[molname] = len(molmoleculedict)
                #         key = str(molmoleculedict[molname] + 1)
                #         if (key in molmolhist):
                #             molmolhist[key] += 1
                #         else:
                #             molmolhist[key] = 1



                elif re.match('End frame', curr):
                    # print(framenum)
                    # for i in molhist:
                    #     molhistperframe[str(framenum) + ',' + i] = molhist[i]
                    newmolhist = molPresent(bond_prev, config_prev, newmoldict)
                    for i in newmolhist:
                        newmolhistperframe[str(framenum) + ',' + str(i)] = int(newmolhist[i])
                    # DumpMolhistperframe(newmolhistperframe, 'newmolhistperframe_' + filetag + '.dat',1)
                    # newmolhistperframe = {}
                    if framenum<2000:
                        for i in config_prev:
                            configperframe[str(framenum) + ',' + str(i)] = config_prev[i]
                        for i in bond_prev:
                            string = []
                            config = []
                            for j in range(len(bond_prev[i])):
                                # string.append(bond_prev[i][j])
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
                    #bond_dist = countBondDist(bond_prev, config_prev, bond_dist, reaction_creat)

                    # for j in bond_prev.keys():
                    #     for k in bond_prev[j]:
                    #         if config_prev[int(k)] > config_prev[int(j)]:
                    #             if str(framenum) + ',(' + str(config_prev[int(j)]) + ' (' + str(config_prev[int(k)]) +' ))' in count_of_bond.keys():
                    #                 count_of_bond[str(framenum) + ',(' + str(config_prev[int(j)]) + ' (' + str(config_prev[int(k)]) +' ))'] += 0.5
                    #             else:
                    #                 count_of_bond[str(framenum) + ',(' + str(config_prev[int(j)]) + ' (' + str(config_prev[int(k)]) +' ))'] = 0.5
                    #         else:
                    #             if str(framenum) + ',(' + str(config_prev[int(k)]) + ' (' + str(config_prev[int(j)]) +' ))' in count_of_bond.keys():
                    #                 count_of_bond[str(framenum) + ',(' + str(config_prev[int(k)]) + ' (' + str(config_prev[int(j)]) +' ))'] += 0.5
                    #             else:
                    #                 count_of_bond[str(framenum) + ',(' + str(config_prev[int(k)]) + ' (' + str(config_prev[int(j)]) +' ))'] = 0.5



                    #REMOVE PART ON MOLHIST, CALCULATE NEWMOLHISTPERFRAME, CALCULATE BOND_DIST ONLY FOR ATOMS WHICH ARE REACTING TOGETHER, COUNT THE NUMBER OF
                    # REACTIONS THAT COULD HAVE HAPPENED FOR DIFFERENT DISTANCES (REACTPERDIST EQUIVALENT)



                    # if framenum<6000:
                    #     string = []
                    #     config = []
                    #     for i in bond_prev:
                    #         string.append(bond_prev[i][j])
                    #         config.append([config_prev[bond_prev[i][j]],bond_prev[i][j]])
                    #         config.sort()
                    #         for j in range(len(bond_prev[i])):
                    #             string.append(config[j][1])
                    #         bond[str(framenum) + ',' + i] = string


                # elif re.match('The box volume', curr):
                #     continue
                #
                # elif re.match('The number of molecules found', curr):
                #     splitline = curr.split(' = ')
                #     molsperframe.append(int(splitline[1]))
                #
                # elif re.match('Beginning molecule', curr):
                #     continue
                #
                # elif re.match('Name', curr):
                #     splitline = curr.split(': ')
                #     splitline2 = splitline[1].split(' \n')
                #     molname = splitline2[0]
                #     if molname not in moleculedict:
                #         moleculedict[molname] = len(moleculedict)
                #
                #         # record number of C's, H's, and O's
                #         # atoms types are currently hard-coded instead of generalized
                #         numC = 0
                #         numH = 0
                #         molnamesplit = molname.split()
                #         j = 0   # i works, as iterator is set, but renamed to avoid confusion
                #         while j < len(molnamesplit) and re.match('(C|H)\d', molnamesplit[j]):
                #             molatoms = re.match('(C|H)(\d*)', molnamesplit[j]).groups()
                #             if molatoms[0] == 'C':
                #                 numC += int(molatoms[1])
                #             elif molatoms[0] == 'H':
                #                 numH += int(molatoms[1])
                #
                #             j += 1
                #         moleculeatomdict.append([numC, numH])
                #
                #         # record number of C-C, C-H, and H-H bonds
                #         # bond types are currently hard-coded instead of generalized
                #         tempbondsdict = {'numCC':'0', 'numCH':'0', 'numHH':'0'}
                #         while j < len(molnamesplit):
                #             bondsplit = molnamesplit[j].split('(')
                #             if bondsplit[1] == 'C-C)':
                #                 tempbondsdict['numCC'] = bondsplit[0]
                #             if bondsplit[1] == 'H-C)':
                #                 tempbondsdict['numCH'] = bondsplit[0]
                #             if bondsplit[1] == 'H-H)':
                #                 tempbondsdict['numHH'] = bondsplit[0]
                #             j += 1
                #         value = tempbondsdict['numCC'] + ' ' + tempbondsdict['numCH'] + ' ' + tempbondsdict['numHH']
                #         bondsdict[molname] = value
                #
                #     try:
                #         key = str(framenum + 1) + ',' + str(moleculedict[molname] + 1)
                #         if(key in molhistperframe):
                #             molhistperframe[key] += 1
                #         else:
                #             molhistperframe[key] = 1
                #     except:
                #         print('filecount: {}'.format(filecount))
                #         print('framenum: {}'.format(framenum))
                #         print('moldictnum: {}'.format(moleculedict[molname]))
                #         raise
                #
                # elif re.match('Atom list', curr):
                #     continue    # do nothing with this for now
                #
                # else:
                #     continue    # do nothing

        # for i in range(framenum):
        #     for j in newmoldict.keys():
        #         if j in moleculedict.keys():
        #             if str(i+1) + ',' + str(moleculedict[j]+1) in molhistperframe.keys():
        #                 newmolhistperframe[str(i+1) + ',' + str(newmoldict[j]+1)] = molhistperframe[str(i+1) + ',' + str(moleculedict[j]+1)]
        #             else:
        #                 newmolhistperframe[str(i + 1) + ',' + str(newmoldict[j] + 1)] = 0
        #         else:
        #             # count = 0
        #             # for k in range(num_atoms):
        #             #     count += checkBondedAtom(j, 0, [[k+1]], 0, bond, configperframe, i+1)
        #             if str(i + 1) + ',' + str(j) in count_of_bond.keys():
        #                 newmolhistperframe[str(i + 1) + ',' + str(newmoldict[j] + 1)] = int(count_of_bond[str(i + 1) + ',' + str(j)])
        #             else:
        #                 newmolhistperframe[str(i + 1) + ',' + str(newmoldict[j] + 1)] = 0





        nummols = len(moleculedict)
        # print('Current number of molecule species observed = {}, filecount = {}'.format(nummols, filecount))
        # DumpMolsperframe(molsperframe, 'molsperframe_' + filetag + '.dat')
        #DumpMolhistperframe(molhistperframe, 'molhistperframe_' + filetag + '.dat',0)
        # DumpMolhistperframe(molmolhistperframe, 'molmolhistperframe_' + filetag + '.dat',0)
        DumpMolhistperframe(newmolhistperframe, 'newmolhistperframe_' + filetag + '.dat',0)
        # DumpBondperframe(count_of_bond, 'bondperframe_' + filetag + '.dat',0)
        # DumpMolhistperframe(configperframe, 'configperframe_' + filetag + '.dat',0)
        # Note molhistperframe will have different number of columns over files, but it's consistent;
        # we can just zero-pad later in MATLAB
    print(bond_dist)
    print("Total number of molecule species observed: {}".format(nummols))
    #DumpMoleculedict(moleculedict, 'moleculedict_all.txt')
    DumpMoleculedict(newmoldict, 'newmoldict.txt')
    DumpMoleculedict(molmoleculedict, 'molmoleculedict_all.txt')
    # DumpMoleculeatomdict(moleculeatomdict, 'moleculeatomdict_all.dat')
    #bonded_reaction = unique(bonded_reaction)


    numreacts = len(newreactionbasisdict)
    print('Current number of distinct reactions observed = {}, filecount = {}'.format(numreacts, filecount))
    #DumpReactlistperframe(reactlistperframe, 'reactperframe_' + filetag + '.dat', numreacts)
    DumpReactlistperframe(newreactlistperframe, 'newreactperframe_' + filetag + '.dat', numreacts)
    #DumpDistanceReactionRate(bonded_reaction, bond_dist, 'distancereactionrate_' + filetag + '.dat', numreacts)

    print('Total number of distinct reactions observed: {}'.format(numreacts))
    #DumpReactionbasisdict(reactionbasisdict, 'reactdict_all.txt')
    DumpReactionbasisdict(newreactionbasisdict, 'newreactdict_all.txt')
    #DumpReactionbasis(reactionbasis, 'reactbasis_all.dat')
    #DumpReactantsbasis(reactantsbasis, 'reactantsbasis_all.dat')
    DumpReactionbasis(newreactionbasis, 'newreactbasis_all.dat')
    DumpReactantsbasis(newreactantsbasis, 'newreactantsbasis_all.dat')
    #DumpNumAtomInvolved(num_atoms_involved, 'num_atoms_involved_all.dat')

    # return [moleculedict, newmoldict, bonded_reaction, bond_dist]
    return




def DumpMolsperframe(molsperframe, dumpmolsfile):
    f = open(dumpmolsfile, 'w')
    for i in range(len(molsperframe)):
        f.write("%s\n" % molsperframe[i])
    f.close()
    print('Done with molecules per frame.')
    return

def DumpMolhistperframe(molhistperframe, dumpmolhistfile,tag):
    if tag ==0:
        f = open(dumpmolhistfile, 'w')
        for i in molhistperframe.keys():
            if molhistperframe[i]!=0:
                [numrow, numcol] = i.split(',')
                f.write("%s %s %s\n" % (numrow, numcol, molhistperframe[i]))
        f.close()
        # print('Done with moleculehist per frame.')
        return
    else:
        f = open(dumpmolhistfile, 'a')
        for i in molhistperframe.keys():
            if molhistperframe[i] != 0:
                [numrow, numcol] = i.split(',')
                f.write("%s %s %s\n" % (numrow, numcol, molhistperframe[i]))
        f.close()
        # print('Done with moleculehist per frame.')
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
    # f2 = open(dumpbondsfile, 'w')
    for key,value in sorted(moleculedict.items(), key = operator.itemgetter(1)):
        f1.write("%s - %s  \n" % (value + 1, key))      # 1-index the molecule labels
        # f2.write("%s\n" % bondsdict[key])   # write corresponding bonds labels
    f1.close()
    # f2.close()
    print('Done with molecule and bonds dictionary.')
    return



def DumpMoleculeatomdict(moleculeatomdict, dumpatomdictfile):
    f = open(dumpatomdictfile, 'w')
    for i in range(len(moleculeatomdict)):
        f.write("%s,%s\n" % (moleculeatomdict[i][0], moleculeatomdict[i][1]))
    f.close()
    print('Done with molecule atom dictionary.')
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

def checkBondedAtom(name, layer, table, index, bond, configperframe, frame):
    l = layer
    t = table
    r = 0
    if index +1 >= len(name):
        return 1
    char = name[index]
    if char == ')':
        l -= 1
        del t[-1]
        r += checkBondedAtom(name, l, t, index + 1, bond, configperframe, frame)
    if char == '(':
        for i in range(len(t[l])):
            if str(configperframe[str(frame) + ',' + str(t[l][i])]) == name[index + 1:index + 5]:
                tt = list(t)
                tt.append(bond[str(frame) + ',' + str(t[l][i])])
                tt[l] = t[l][i+1:len(t[l])]
                r += checkBondedAtom(name, l + 1, tt, index + 6, bond, configperframe,
                              frame)
    return r

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

def countBondDist(bond, config, bond_dist, reaction_creat):
    bond_dist_new = bond_dist
    visited = [0]*len(bond)
    for i in range(len(bond)):
        if visited[i] != 1:
            atoms_in_mol = getMolecule(i+1, bond, [])
            bond_dist_new = countBondDistMol(atoms_in_mol, bond, config, bond_dist_new, reaction_creat)
            for j in atoms_in_mol:
                visited[j-1] = 1
    return bond_dist_new

def getMolecule(index, bond, atoms_in_mol):
    if index in atoms_in_mol:
        return atoms_in_mol
    else:
        atoms_in_mol.append(index)
        for i in bond[str(index)]:
            atoms_in_mol = getMolecule(i, bond, atoms_in_mol)
    return atoms_in_mol

def countBondDistMol(atoms, bond, config, bond_dist_new, reaction_creat):
    for i in range(len(atoms)-1):
        for j in range(i+1,len(atoms)):
            c1 = int(config[atoms[i]])
            c2 = int(config[atoms[j]])
            k = are_bonded(atoms[i], atoms[j], [], bond)
            if c1>c2:
                for l in range(len(reaction_creat)):
                    if str(c2) == reaction_creat[l][0] and str(c1) == reaction_creat[l][1]:
                        if k in bond_dist_new.keys():
                            bond_dist_new[k] +=1
                        else:
                            bond_dist_new[k] = 1
            else:
                for l in range(len(reaction_creat)):
                    if str(c2) == reaction_creat[l][1] and str(c1) == reaction_creat[l][0]:
                        if k in bond_dist_new.keys():
                            bond_dist_new[k] +=1
                        else:
                            bond_dist_new[k] = 1
    return bond_dist_new

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


def DumpDistanceReactionRate(bonded_reaction, reactperdist, dumpreactperframe, numreacts):
    f = open(dumpreactperframe, 'w')
    for i in reactperdist.keys():
        if int(i) in bonded_reaction[0]:
            for j in range(len(bonded_reaction[1])):
                if int(i) == bonded_reaction[0][j]:
                    index = j
            f.write("%s %s %s %s\n" % (
            i, reactperdist[i], bonded_reaction[1][index], bonded_reaction[1][j] / reactperdist[i]))
        else:
            f.write(
                "%s %s %s %s\n" % (i, reactperdist[i], 0, 0))
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


def DumpNumAtomInvolved(reactionbasisdict, dumpdict):
    f = open(dumpdict, 'w')
    for key in reactionbasisdict.keys():
        f.write("%s %s\n" % (key + 1, reactionbasisdict[key]))
    f.close()
    print('Done with reaction basis dictionary.')
    return


def getBondedAtom(Left, Left_Bond):
    length = [len(Left), len(Left_Bond)]
    unique = set(length)
    if (len(unique) != 1):
        print(length)
        print('Error: Not same length for Left, Right, Left_Bond or Right_Bond')
        return
    reactant = []
    store = {}
    for i in range(len(Left)):
        split = Left_Bond[i].split(' ')
        for j in range(len(split) - 1):
            split[j] = int(split[j])
        if split[0] == 0:
            reactant.append(Left[i])
        else:
            store[i + 1] = []
            for k in range(len(split) - 1):
                store[i + 1].append(split[k])
    visited = {}
    for k in store.keys():
        visited[k] = 0
    for k in store.keys():
        string = ''
        string = StringBondedAtom(string, k, store, visited, Left)
        if string != '':
            reactant.append(string)
    return reactant


def StringBondedAtom(string, m, store, visited, atom_list):
    if visited[m] == 1:
        return string
    string += '(' + str(atom_list[m - 1])
    visited[m] = 1
    for i in range(len(store[m])):
        string = StringBondedAtom(string, store[m][i], store, visited, atom_list)
    string += ')'
    return string




