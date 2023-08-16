import numpy as np
from ase.build import bulk
from ase.io import write,read
import random
from ase import Atoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.neighborlist import neighbor_list
import time
from datetime import date, time, datetime

def number_det(atom,elements,crystal):
    '''

    '''
    num_types = len(elements)
    #for ionic structure
    
    if (crystal == 'rocksalt'):
        num_atom = int (len(atom)/2) #as half of the atoms are anion
    elif (crystal == 'fcc' or crystal == 'bcc' or crystal == 'hcp'):
        num_atom = len(atom)
    else:
        raise Exception("Crystal structure not implemented yet!")

    num_elem = num_atom/num_types

    #checking if total number of elements id divisible by type of elements
    check_int = isinstance(num_elem,int)
    num_each = []
    if(check_int == False):
        near_int = round(num_elem) 
        for a in range(num_types):
            if(a == (num_types-1)):
                extra_num = near_int + (num_atom - (num_types * near_int))
                # By above, we make sure that whole number of atoms are for each species
                # which leads to an extra atom for last index 
                num_each.append(extra_num)
            else:
                num_each.append(near_int)
    if(check_int == True):
        for a in range(num_types):
            num_each.append(num_elem)

    return num_each
def element_assign(elements,num_each):
    '''
    This module calculates the list of atoms. 
    Such list contains the the name of chemical
    species in a sequence. 
    Input:
    elements: list containing chemical species
    num_each: list containing number of each
    chemical species.
    Returns:
    elem_list : List of atoms in a sequence, such 
    that, same elements are together. 
    '''
    elem_list = []
    for a,b in enumerate(num_each):
        for c in range(b):
            elem_list.append(elements[a])
    
    return elem_list
def random_gen(metal_pos,atom):
    '''
    This function simply randomises the metal positions
    which are basically all positions for alloys, while
    these are simply cation positions for ionic solids.
    '''
    #num_list = []
    #for a in range(len(atom)):
    #    num_list.append(a)
    
    random.shuffle(metal_pos)
    positions = atom.get_positions()
    rand_pos = []
    for a in range(len(metal_pos)):
        rand_pos.append(positions[metal_pos[a]])

    return rand_pos

def elem_dict(elements):
    '''

    '''
    elem_id = {}
    for a,b in enumerate(elements):
        elem_id[b] = a + 1
    return elem_id
def cutoff_det(d,crystal,fluc_fact):
    ctoff = 0.
    if(crystal=="fcc"):
        ctoff = (d/(np.sqrt(2))) + fluc_fact			
        return ctoff
    elif(crystal=="bcc"):
        ctoff = ((np.sqrt(3)*d)/(2.)) + fluc_fact  
        return ctoff
    elif(crystal=="hcp"):
        ctoff = d + fluc_fact							
        return ctoff
    elif (crystal == 'rocksalt'):
        #1NN for rocksalt is d/2, while 2NN is d/sqrt(2)
        ctoff = d/1.414 + fluc_fact
        return ctoff
    else:
        raise Exception('Unknown input to cutoff det function')
def maxnn_det(crystal):
    if(crystal=="bcc"):
        nn = 8
        return nn
    elif(crystal=="fcc" or crystal=="hcp"):
        nn = 12
        return nn
    elif (crystal == 'rocksalt'):
        nn = 18
        return nn
    else:
        raise Exception('Unknown input to maxnn_det function')

def num_det(elements):
    num=int((len(elements)*(len(elements)+1))/2)
    return num

def pairs_gen(elements,num_pairs):
    pairs=[]
    for a in range(len(elements)):
        for b in range(a,len(elements)):
            pairs.append('{}-{}'.format(elements[a],elements[b]))
    return pairs

def bond_pair_count(pairs,i,j,elem_list):
    counter_val1=np.zeros(len(pairs)) #array to store the no of bonds.
    counter_val2=np.zeros(len(pairs))
    counter_val3=np.zeros(len(pairs))
    #generating the bond_pair list
    bond_pair=[]
    #print ('{},{},{}'.format(len(i),len(j),len(elem_list))) #TODO remove
    for a in range(len(i)):
        #print ('{},{}'.format(i[a],j[a])) #TODO remove this
        bond_pair.append('{}-{}'.format(elem_list[i[a]],elem_list[j[a]]))
    for b,c in enumerate(pairs):
        counter_val1[b]=bond_pair.count('{}-{}'.format(c.split('-')[0],c.split('-')[1]))
        counter_val2[b]=bond_pair.count('{}-{}'.format(c.split('-')[1],c.split('-')[0]))
        if (c.split('-')[0] != c.split('-')[1]): #unlike bond
            counter_val3[b]=counter_val1[b]+counter_val2[b]
        else:
            counter_val3[b]=counter_val1[b] 
    #end9=time.time()
    #print ('time-taken9,counter_val={},{}'.format(end9-start9,counter_val3))
    #counter_val=np.zeros(len(pairs)) #array to store the no of bonds. #TODo remove after test
    #for a in range(len(i)):
    #    bond_pair1='{}'.format(elem_list[i[a]])
    #    bond_pair2='{}'.format(elem_list[j[a]])
    
    #    for b,c in enumerate(pairs):
            #print ('{}-{}'.format(bond_pair1,bond_pair2),c)
    #        if(('{}-{}'.format(bond_pair1,bond_pair2) == c) or \
    #                ('{}-{}'.format(bond_pair2,bond_pair1)) == c):
    #            counter_val[b]=counter_val[b]+1
    #print ('counter-value after={}'.format(counter_val)) #TODO remove this after test
    return counter_val3

def req_bond_val(sro,pairs,pref_pair,elements,i,atoms,crystal):
    req_val=np.zeros(len(pairs))
    #finding the index in pairs list for pre_pair
    count=0
    indx=-1
    for item in pairs:
        temp=[]
        temp.append(item.split('-'))
        #print (len(temp),temp,pref_pair[0],pref_pair[1])
        #print (count,temp[0],pref_pair[0],temp[1],pref_pair[1])
        if((temp[0][0] == pref_pair[0]) and (temp[0][1] == pref_pair[1])):
            indx=count
            break
        elif ((temp[0][1] == pref_pair[0]) and (temp[0][0] == pref_pair[1])):
            indx=count
            break
        else:
            count+=1
def indx_find(num_each,num,elements):
    temp=0
    val=-1
    for a,b in enumerate(elements):
        if(num < temp+num_each[a]):
            val=a
            break
        else:
            temp+=num_each[a]
    return val
def order_para_cal(pairs,count_bonds,num_like,num_unlike,total_num_bonds):
    #determination of number of like bonds in the scenarion of complete disorder(num_bond_disorder).
    num_bond_disorder=(total_num_bonds)/((2*num_unlike)+num_like)
    unlike_counter=0.
    like_counter=0.
    for num,item in enumerate(pairs):
        temp1=item.split('-')[0]
        temp2=item.split('-')[1]
        if(temp1 != temp2):
            #TODO here num_like = num_unlike bonds are being assumed at complete disorder
            unlike_counter=((1.-(count_bonds[num]/(2*num_bond_disorder)))/num_unlike)+unlike_counter
        elif(temp1 == temp2):
            like_counter=(((count_bonds[num]/num_bond_disorder)-1.)/num_like)+like_counter
        else:
            print('Issue with OP calculation')
            pass
    delta=unlike_counter+like_counter
    return delta,num_bond_disorder
        
def pair_index(pairs,at1,at2):
    for num,item in enumerate(pairs):
        temp1=item.split('-')[0]
        temp2=item.split('-')[1]
        if(temp1 == at1 and temp2 == at2):
            break
        elif(temp1 == at2 and temp2 == at1):
            break
        else:
            pass
    return num
def nn_dict_det(i,j):
    '''
    Parameter:
        i: Atom identity
        j: Nearest-neighbour (NN) of ith atom

    Returns:
        nn_dict: Dictionary with key (i) and values (NN atoms of i)
    '''
    nn_dict={}
    item_ear=-1
    for num,item in enumerate(i):
        if(item == item_ear):
            nnlist.append(j[num])
        elif(item != item_ear and num > 0 and num < len(i)):
            nn_dict[item_ear]=nnlist
            nnlist=[]
            nnlist.append(j[num])
            item_ear=item
        elif(item != item_ear and num == 0):
            nnlist=[]
            nnlist.append(j[num])
            item_ear=item
    nn_dict[item_ear]=nnlist
    return nn_dict

def fl_check(forbidden_list,a):
    check=False
    if(len(forbidden_list) == 0):
        pass
    else:
        #check=[True for item in forbidden_list if item == a]
        for item in forbidden_list:
            if(item == a):
                check=True
                break
            else:
                pass
    return check
def count_proxy(proxy_ele,list_element):
    counter=0
    counter=sum(1 for item in list_element if item == proxy_ele)
    return counter
def forbidden_list_count(forbidden_list):
    counter=0
    counter=sum(1 for item in forbidden_list if item != -1)
    return counter
def ele_indx(elements,val):
    '''

    '''
    indx = -1
    for num,item in enumerate(elements):
        if (item == val):
            indx = num
            break
        else:
            pass
    return indx
def start_end(elements,val,num_each):
    '''

    '''
    ele_id = ele_indx(elements,val)
    start = 0
    end = 0
    for a in range(len(elements)):
        if (a < ele_id):
            start = start + num_each[a]
            end = end + num_each[a]
        elif (a == ele_id):
            end = end + num_each[a]
        else:
            pass
    return start,end
def lat_para_det(ele_rad,an_rad,crystal):
    '''
    
    '''
    if (crystal == 'fcc'):
        rmean = 0.
        for val in ele_rad:
            rmean = rmean + val
        rmean = rmean/len(ele_rad)
        
        lat_para = 2.*1.414*rmean
        return lat_para
    elif (crystal == 'bcc'):
        rmean = 0.
        for val in ele_rad:
            rmean = rmean + val;
        rmean = rmean/len(ele_rad)
        
        lat_para = (4*rmean)/1.717
        return lat_para
    elif (crystal == 'rocksalt'):
        rmean = 0.
        for val in ele_rad:
            rmean = rmean + val;
        rmean = rmean/len(ele_rad)
        lat_para = (2*(an_rad + rmean))/1.717
        return lat_para
    elif (crystal == 'hcp'):
        rmean = 0.
        for val in ele_rad:
            rmean = rmean + val
        rmean = rmean/len(ele_rad)
        lat_para = 2*rmean
        return lat_para
    else:
        raise Exception('Crystal structure not implemented')

def log_write(cal_mode,pairs,count_bonds,delta,end):
    '''
    
    '''
    if (cal_mode == 0):
        pass #we are not writing log file for the disordered structure generation.
    if (cal_mode == 1):
        with open('OPERA.log','a') as out:
            out.write(' '.join(pairs))
            out.write('\t')
            for val in count_bonds:
                out.write('{}'.format(str(val)))
                out.write(' ')
            out.write('\t')
            out.write(str('{:.4f}'.format(delta)))
            if (end == True):
                out.write('\n')
                out.write('log file written at:{}'.format(datetime.now()))
                out.write('\n')
            elif (end == False):
                out.write('\n')
