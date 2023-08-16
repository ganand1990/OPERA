import func_all
import numpy as np
from ase.io import write,read
from ase import Atoms
from ase.neighborlist import neighbor_list
import random_gen
import max_swap
import os,glob
import time
import bib
import warnings
def calculation_mode():
    cal_mode = 1#calculation mode; 0 for random structure gen and 1 for reading the random file and generating SRO confs
    return cal_mode
def input_data():
    #elements = ['Sc','Ti','Zr','Hf'] #this order is dependent upon the ordering the potential file
    #ele_rad = [1.62,1.47,1.60,1.59]
    #elements = ['Mo','Nb','Ta','W']
    #ele_rad = [1.39,1.46,1.46,1.39]
    elements = ['Co','Cr','Ni']#'Fe','Mn','Ni']#,'Ni'] #Cr removed
    ele_rad = [1.25,1.28,1.24]#,1.24] #Cr:1.28 Ang, Fe: 1.26, Mn: 1.27
    anion = ['O'] #add any anion (not same as in element list) as proxy here.
    an_rad = 0. #Anion radius; Set it zero for non-ionic structure. 
    crystal='fcc'
    Tmax = 5 #parameter for max temerature in simulated annealing
    Tnum = 1000 #parameter for number of steps for T -> 0 in SA
    #order in which, atoms are named in pref_pair does not matter.
    fluc=0.1
    inp_file='inp.cfg' #input file for random structure generation
    random_file = '0.0_random.xyz'
    proxy_ele='Ca'
    pref_pair_req=['Co','Co'] #it can be in the arbitrary order
    num_delta=5 #number of delta para from 0 to max delta possible
    return elements,anion,ele_rad,an_rad,crystal,fluc,inp_file,random_file,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum

def make_data(calc_mode):
    elements,anion,ele_rad,an_rad,crystal,fluc,inp_file,random_file,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum=input_data()
    #Guessing the lattice parameter with Vegard's Law
    lat_param = func_all.lat_para_det(ele_rad,an_rad,crystal)
    if ( calc_mode == 0 ):    
        at = read('{}'.format(inp_file))
    elif ( calc_mode == 1 ):
        at = read('{}'.format(random_file))
    else:
        raise Exception('Unidentified calc_mode value')

    #TEST for ionic structure, whether the inp.cfg has same anion as desired.
    if (crystal == 'rocksalt'):
        anion_pres = (at.symbols == anion)
        anion_val = sum(1 for val in anion_pres if val == True)
        if (anion_val < 1):
            raise Exception('anion in inp.cfg is not same as anion as desired.')
    elif (crystal == 'fcc' or crystal == 'bcc' or crystal == 'hcp'):
        #if (len(anion) != 0):
        #    raise Exception('anion list should be empty')
        pass
    else:
        raise Exception('crystal structure not implemented')
    #Note that in case of alloys unchangable sublattice is not there
    #num_each list determines the presence or absence of unchangable 
    #sublattice. For e.g., for oxides the sum of num_each is number of 
    #cations, while for alloys it is simply number of atoms.

    #determining number of each chemical species
    num_each = func_all.number_det(at,elements,crystal)
            
    #list of atoms in a sequence
    elem_list = func_all.element_assign(elements,num_each)
    #position of anions and cations are stored separately
    #for alloys, anion array would be empty.
    anion_pres = (at.symbols == anion)
    anion_pos = [a for a,b in enumerate(anion_pres) if b == True]
    total_pos = [a for a in range(len(at))]
    anion_pos_set = set(anion_pos)
    total_pos_set = set(total_pos)
    metal_pos = list(total_pos_set.symmetric_difference(anion_pos_set))
    all_pos = at.get_positions()
    anion_coord = []
    for item in anion_pos:
        anion_coord.append(all_pos[item])
    
    #counting number of pairs of like atoms
    num_pairs=func_all.num_det(elements)
    num_like=len(elements)
    num_unlike=num_pairs-num_like
    pairs=[]
    #generating bonding pairs
    pairs=func_all.pairs_gen(elements,num_pairs)
    cutoff=func_all.cutoff_det(lat_param,crystal,fluc)

    if ( calc_mode == 0 ):
        rand_pos = func_all.random_gen(metal_pos,at)
    elif ( calc_mode == 1):
        rand_pos = [] #position of metals
        for item in metal_pos:
            rand_pos.append(all_pos[item])
    else:
        raise Exception('Unidentified calc_mode value')


    if (calc_mode == 0 ):
        atoms = Atoms(cell=at.get_cell())
        #putting metal at their positions
        for a in range(len(rand_pos)):
            atoms.extend(Atoms('{}'.format(elem_list[a]),positions=[(rand_pos[a][0],rand_pos[a][1],rand_pos[a][2])]))
        for a in range(len(anion_coord)):
            atoms.extend(Atoms('{}'.format(anion[0]),positions=[(anion_coord[a][0],anion_coord[a][1],anion_coord[a][2])]))
        atoms.set_pbc(111)
    elif ( calc_mode == 1 ): 
        atoms = at
        num_elements = len(set(atoms.get_chemical_symbols())) 
        '''
        implementing a condition that if the random structure file
        for different number of elements is being read
        '''
        if (num_elements != len(elements)):
            print ('num_elements,len(elements)={},{}'.format(num_elements,len(elements)))
            raise Exception ('Random structure file has diffent no. of elements wrt element list of the present calculation')
        else:
            pass
    else:
        raise Exception('Unindentified calc_mode value')
    
    
    #cutoff=cutoff_det(latparam,crystal,fluc)
    i,j = neighbor_list('ij',atoms,cutoff)
    #modify i and j to include only metals and remove anions
    i_mod = []
    j_mod = []
    for a,b in zip(i,j):
        if ( a < len(rand_pos) and b < len(rand_pos)):
            i_mod.append(a)
            j_mod.append(b)
    cn=int (len(i)/len(at))
    #print (cn)
    #TEST to ensure the right CN
    cnmax=func_all.maxnn_det(crystal)
    if(cnmax < cn):
        raise Exception("Decrease fluc_fact to make CN right!\n")
    if(cnmax > cn):
        raise Exception("Increase the fluc_fact to make CN right!\n")
    
    return elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req\
            ,num_delta,Tmax,Tnum,i_mod,j_mod,elem_list,rand_pos,anion_coord,anion
def div_equal(num_delta,num_swap_pos):
    '''

    '''
    list_num_swap=np.full(num_delta,-1,dtype=int)
    val=int (num_swap_pos/num_delta)
    count=0
    ear_val=0
    for num in range(num_delta):
        if(num == num_delta-1):
            list_num_swap[num]=val + (num_swap_pos-(count+val)) + ear_val
        else:
            list_num_swap[num]=val + ear_val
            count+=val
            ear_val=list_num_swap[num]
    return list_num_swap
def atoms_gen(num_swap,atoms,swap_dict):
    '''
        Parameters:
            num_swap:
            atoms:
            swap_dict:
        Returns:
            at:
    '''
    counter=0
    elem_list=atoms.get_chemical_symbols()
    pos=atoms.get_positions()
    at=Atoms(cell=atoms.get_cell())
    for item in swap_dict:
        if (counter < num_swap):
            first=np.full(3,0,dtype=float)
            second=np.full(3,0,dtype=float)
            first=[a for a in pos[item]]
            second=[a for a in pos[swap_dict[item]]]
            #print ('swap_check={},{},{},{}'.format(pos[item],pos[swap_dict[item]],first,second))
            #SWAP
            #print ('temp1,temp2={},{}'.format(first,second))
            #pos[item]=[a for a in second]
            #pos[swap_dict[item]]=[b for b in first]
            for num,a in enumerate(second):
                pos[item][num]=a
            for num,b in enumerate(first):
                pos[swap_dict[item]][num]=b
            #print ('after temp1,temp2={},{}'.format(first,second))
            #pos[item]=temp2
            #pos[swap_dict[item]]=temp1
            #print ('after swap_check={},{},{},{}'.format(pos[item],pos[swap_dict[item]],first,second))
            counter+=1
    for a in range(len(atoms)):
        at.extend(Atoms('{}'.format(elem_list[a]),positions=[(pos[a][0],pos[a][1],pos[a][2])]))
    at.set_pbc((True,True,True))
    return at
def sro_count_bond(atoms,cutoff,pairs):
    '''
        Parameters:
            atoms:
            cutoff:
            pairs:
        Returns:
            count_bond:
    '''
    #determination of nn and nn_dict,bond count 
    i,j=neighbor_list('ij',atoms,cutoff)
    nn_dict=func_all.nn_dict_det(i,j)
    elem_list=atoms.get_chemical_symbols()
    positions=atoms.get_positions()
    count_bond=func_all.bond_pair_count(pairs,i,j,elem_list)
    return count_bond
def bond_count_trend(pref_pair_req,pairs,elements,num_each,nn_dict,count_bond,num_delta,atoms,cutoff):
    '''
        Parameters::
            pref_pair_req:
            pairs:
            elements:
            num_each:
            nn_dict:
            count_bond:
            num_delta:
            atoms:
            cutoff:

        Returns:
            bond_trend:
            swap_dict:
            num_swap_pos:

    '''
    num_swap_pos,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2=max_swap.swap_pos(pref_pair_req,\
                pairs,elements,num_each,nn_dict,count_bond)
    #print ('num_swap_pos,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2={},{},{},{},{},{}'.format(\
    #        num_swap_pos,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2))
    list_num_swap=np.full(num_delta,-1,dtype=int)
    list_num_swap=div_equal(num_delta,num_swap_pos)
    #print ('list_num_swap={}'.format( list_num_swap))
    if (num_swap_pos < num_delta):
        #this condition can arise if the num_delta is greater than max swap possible.
        bond_trend = False
        warnings.warn('If this persists, decrease the num_delta value or increase the number of atoms in the system')
    else:
        bond_trend=True
        '''
        A check is being introduced that swap leads to the increase in the desired bond and 
        decrease in the undesired bonds.
        '''
        for num,item in enumerate(list_num_swap):
            at_ear=atoms
            at=atoms_gen(item,at_ear,swap_dict)
            count_bonds_re=sro_count_bond(at,cutoff,pairs)
            #TODO There is chance of modification below to tailor the increase or decrease 
            #of the desired and undesired bonds, respectively.
            if ((count_bonds_re[indx_red1] > count_bond[indx_red1]) or \
                (count_bonds_re[indx_red2] > count_bond[indx_red2]) or \
                (count_bonds_re[indx_inc1] < count_bond[indx_inc1]) or \
                (count_bonds_re[indx_inc2] < count_bond[indx_inc2])):
                bond_trend=False
                break
            else:
                count_bond=count_bonds_re #TODO seems redundant!
    return bond_trend,swap_dict,num_swap_pos    

def main():
    cal_mode = calculation_mode()
    if( cal_mode == 0 ):
        start=time.time()
        elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,\
                pref_pair_req,num_delta,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord\
                ,anion=make_data(cal_mode)
        count_bonds,atom_left,delta,ideal_num=random_gen.random_struc_gen(proxy_ele,\
                elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,Tmax,Tnum,\
                i,j,elem_list,rand_pos,anion_coord,anion)
        print ('pairs,count_bonds,delta={},{},{}'.format(pairs,count_bonds,delta))
        end=time.time()
        print ('time-taken={}'.format(end-start))
        #printing citation information
        bib.citation()
    elif ( cal_mode == 1):
        elements,num_each,atoms,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,\
                pref_pair_req,num_delta,Tmax,Tnum,i,j,elem_list,rand_pos,\
                anion_coord,anion=make_data(cal_mode)
        
        
        total_num_bonds=len(i) #total number of bonds
        nn_dict=func_all.nn_dict_det(i,j)
        #Test that if the serial of elements in elem_list and elements ][ should be same.
        pre_ele=elem_list[0]
        set_elem_list=[]
        set_elem_list.append(pre_ele)
        for a,b in enumerate(elem_list):
            if (elem_list[a] != pre_ele):
                set_elem_list.append(elem_list[a])
                pre_ele=elem_list[a]
            else:
                pass
        #print ('set_elem_list={}'.format(set_elem_list))
        if(set_elem_list == elements):
            pass
        else:
            raise Exception('elem_list and elements list dont have same order')
        #positions=atoms.get_positions()
        count_bond=func_all.bond_pair_count(pairs,i,j,elem_list)
        #print('count_bonds={}'.format(count_bond))
        #print('pref_pair_req,pairs,elements,num_each={}{}{}{}'.format(pref_pair_req,pairs,elements,num_each))
        bond_trend=False
        while (bond_trend == False):
            bond_trend,swap_dict,num_swap_pos=bond_count_trend(pref_pair_req,pairs,elements,num_each,\
                    nn_dict,count_bond,num_delta,atoms,cutoff)
            #print ('bond_trend in main={}'.format(bond_trend))
        #generating atoms object with diff delta paramter
        list_num_swap=np.full(num_delta,-1,dtype=int)
        list_num_swap=div_equal(num_delta,num_swap_pos)
        #print (list_num_swap)
        #cleaning the directory

        for filename in glob.glob("*_woutdis.xyz"):
            os.remove(filename)
        for filename in glob.glob("*_woutdis.vasp"):
            os.remove(filename)

        for num,item in enumerate(list_num_swap): 
            at_ear=atoms
            at=atoms_gen(item,at_ear,swap_dict)
            count_bonds_re=sro_count_bond(at,cutoff,pairs)
            delta,num_bond_disorder=func_all.order_para_cal(pairs,count_bonds_re,num_like,num_unlike,total_num_bonds) 
            print (pairs,count_bonds_re,delta)
            func_all.log_write(cal_mode,pairs,count_bonds_re,delta,end=False)
            if(delta < 0.):
                write('neg_{:.4f}_woutdis.xyz'.format(np.abs(delta)),at)
                write('neg_{:.4f}_woutdis.vasp'.format(np.abs(delta)),at)
            else:
                write('pos_{:.4f}_woutdis.xyz'.format(delta),at)
                write('pos_{:.4f}_woutdis.vasp'.format(delta),at)
        #writing the disordered file, since pos value is assigned to 0.0 delta param
        if(delta < 0.):
            write('neg_{:.4f}_woutdis.xyz'.format(0.),atoms)
            write('neg_{:.4f}_woutdis.vasp'.format(0.),atoms)
        else:
            write('pos_{:.4f}_woutdis.xyz'.format(0.),atoms)
            write('pos_{:.4f}_woutdis.vasp'.format(0.),atoms)
        #printing citation information
        bib.citation()
        #printing log file
        func_all.log_write(cal_mode,pairs,count_bonds_re,delta,end=True)
    else:
        raise Exception('cal_mode val not right')
if __name__=="__main__":
    main()
