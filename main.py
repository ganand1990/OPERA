import func_all
import numpy as np
from ase.io import write,read
from ase import Atoms
from ase.neighborlist import neighbor_list
import random_gen
import max_swap
import os,glob
import time
def calculation_mode():
    cal_mode = 1 #calculation mode; 0 for random structure gen and 1 for reading the random file and generating SRO confs
    return cal_mode
def input_data():
    elements = ['Sc','Ti','Zr','Hf'] #this order is dependent upon the ordering the potential file
    ele_rad = [1.62,1.47,1.60,1.59]
    #elements = ['Mo','Nb','Ta','W']
    #ele_rad = [1.39,1.46,1.46,1.39]
    #elements = ['Co','Cr','Ni']
    #ele_rad = [1.246,1.280,1.241]
    anion = ['O'] #add any anion (not same as in element list) as proxy here.
    an_rad = 0. #Anion radius; Set it zero for non-ionic structure. 
    crystal='fcc'
    Tmax = 5 #parameter for max temperature in simulated annealing
    Tnum = 1000 #parameter for number of steps for T -> 0 in SA
    #order in which, atoms are named in pref_pair does not matter.
    fluc=0.1
    inp_file='inp.cfg' #input file for random structure generation
    random_file = '0.0_random.xyz'
    proxy_ele='Ca'
    pref_pair_req=['Sc','Ti'] #it can be in the arbitrary order
    num_delta=5 #number of delta para from 0 to max delta possible
    #cal_mode=1 #calculation mode; 0 for random structure gen and 1 for reading the random file and generating SRO confs
    return elements,anion,ele_rad,an_rad,crystal,fluc,inp_file,random_file,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum

def make_data(calc_mode):
    elements,anion,ele_rad,an_rad,crystal,fluc,inp_file,random_file,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum=input_data()
    #Guessing the lattice parameter with Vegard's Law
    lat_param = func_all.lat_para_det(ele_rad,an_rad,crystal)
    #latparam = 0.
    #for item in latpara_list:
    #    latparam = latparam + item
    #latparam = latparam/(len(latpara_list))
    print (lat_param)
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
    #num_metal = sum(a for a in num_each)
    #num_anion = len(at) - num_metal
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
        print ("Decrease fluc_fact to make CN right!\n")
        exit(1)
    if(cnmax > cn):
        print ("Increase the fluc_fact to make CN right!\n")
        exit(1)
    return elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum,i_mod,j_mod,elem_list,rand_pos,anion_coord,anion
def div_equal(num_delta,num_swap_pos):
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
    #determination of nn and nn_dict,bond count 
    i,j=neighbor_list('ij',atoms,cutoff)
    nn_dict=func_all.nn_dict_det(i,j)
    elem_list=atoms.get_chemical_symbols()
    positions=atoms.get_positions()
    count_bond=func_all.bond_pair_count(pairs,i,j,elem_list)
    return count_bond
def bond_count_trend(pref_pair_req,pairs,elements,num_each,nn_dict,count_bond,num_delta,atoms,cutoff):
    num_swap_pos,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2=max_swap.swap_pos(pref_pair_req,\
                pairs,elements,num_each,nn_dict,count_bond)
    list_num_swap=np.full(num_delta,-1,dtype=int)
    list_num_swap=div_equal(num_delta,num_swap_pos)
    #count_bond_re=counts
    bond_trend=True
    for num,item in enumerate(list_num_swap):
        #putting a condition that, if the list_num_swap has zeros, then bond_trend
        #should be False 
        set_num_swap=set(list_num_swap)
        counter_num_swap=0
        for a in range(len(set_num_swap)):
            counter_num_swap+=1
        if (counter_num_swap < num_delta):
            print('If this persists, decrease the num_delta value or increase the number of atoms in the system')
            #this condition can arise if the num_delta is greater than max swap possible.
            bond_trend=False
            break
        else:
            pass
        at_ear=atoms
        at=atoms_gen(item,at_ear,swap_dict)
    #        at_ear=at
    #    else:
    #        at=atoms_gen(item,at_ear,swap_dict)
    #        at_ear=at
        count_bonds_re=sro_count_bond(at,cutoff,pairs)
        if ((count_bonds_re[indx_red1] > count_bond[indx_red1]) or \
                (count_bonds_re[indx_red2] > count_bond[indx_red2]) or \
                (count_bonds_re[indx_inc1] < count_bond[indx_inc1]) or \
                (count_bonds_re[indx_inc2] < count_bond[indx_inc2])):
            bond_trend=False
            break
        else:
            count_bond=count_bonds_re
    return bond_trend,swap_dict,num_swap_pos    

def main():
    cal_mode = calculation_mode()
    #elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion=make_data()
    if( cal_mode == 0 ):
        start=time.time()
        #random_reach=False
        #while (random_reach == False):
        elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion=make_data(cal_mode)
        #elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req,num_delta,cal_mode,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion=make_data()
        count_bonds,atom_left,delta,ideal_num=random_gen.random_struc_gen(proxy_ele,elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion)
        print ('pairs,count_bonds,delta={},{},{}'.format(pairs,count_bonds,delta))
            #random_reach_count=0
            #for num,item in enumerate(pairs):
            #    if (item.split('-')[0] == item.split('-')[1] and count_bonds[num] == ideal_num):
            #        random_reach_count+=1
            #    elif (item.split('-')[0] != item.split('-')[1] and count_bonds[num] == 2*ideal_num):
            #        random_reach_count+=1
            #    else:
            #        pass
            #if (random_reach_count == len(pairs)):
            #    random_reach=True
            #else:
            #    random_reach=False
        end=time.time()
        print ('time-taken={}'.format(end-start))
    elif ( cal_mode == 1):
        elements,num_each,atoms,cutoff,pairs,num_unlike,num_like,cn,proxy_ele,pref_pair_req,num_delta,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion=make_data(cal_mode)
        
        #reading the atoms object with zero delta paramter
        #atoms=read('0.0_random.xyz')
        #elem_list_all=atoms.get_chemical_symbols()
        #num_metal = sum(a for a in num_each)
        #TEST to ensure that metals are before anion in the elem_list_all
        #generation of elem_list
        #for a,b in enumerate(elem_list_all):
        #    print (a,b)
        #    if ( a < num_metal ):
        #        elem_list.append(b)
        #    else:
        #        pass
        
        #removing the anions from the elem_list
        #determination of nn and nn_dict,bond count 
        #i,j=neighbor_list('ij',atoms,cutoff)
        #i_mod and j_mod need to be calculated, which would contain only metal not anions
        #i_mod = []
        #j_mod = []
        #for a,b in zip(i,j):
        #    if ( a <  num_metal and b < num_metal ):
        #        i_mod.append(a)
        #        j_mod.append(b)

        total_num_bonds=len(i) #total number of bonds
        #print (num_metal,len(elem_list),total_num_bonds)
        
        nn_dict=func_all.nn_dict_det(i,j)
        #print (nn_dict)
        #elem_list=atoms.get_chemical_symbols()
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
        print (set_elem_list)
        if(set_elem_list == elements):
            pass
        else:
            print ('elem_list and elements list dont have same order')
            exit(1)
        #positions=atoms.get_positions()
        count_bond=func_all.bond_pair_count(pairs,i,j,elem_list)
        print(count_bond)
        print(pref_pair_req,pairs,elements,num_each)
        bond_trend=False
        while (bond_trend == False):
            bond_trend,swap_dict,num_swap_pos=bond_count_trend(pref_pair_req,pairs,elements,num_each,\
                    nn_dict,count_bond,num_delta,atoms,cutoff)
            print (bond_trend)
        #1/0
        #num_swap_pos,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2=max_swap.swap_pos(pref_pair_req,\
        #        pairs,elements,num_each,nn_dict,count_bond)
        #print (elem_list)
        #initialising the count_bond_re numpy array
        #count_bond_re=np.full(len(pairs),0,dtype=float)
        #for item in swap_dict:
        #    print (item,swap_dict[item],elem_list[item],elem_list[swap_dict[item]])
        #generating atoms object with diff delta paramter
        list_num_swap=np.full(num_delta,-1,dtype=int)
        list_num_swap=div_equal(num_delta,num_swap_pos)
        print (list_num_swap)
        #cleaning the directory
        #if os.path.exists("*_woutdis.xyz"):
        for filename in glob.glob("*_woutdis.xyz"):
            os.remove(filename)
        for filename in glob.glob("*_woutdis.vasp"):
            os.remove(filename)

        for num,item in enumerate(list_num_swap):
        #    if(num==0):
            at_ear=atoms
            at=atoms_gen(item,at_ear,swap_dict)
        #        at_ear=at
        #    else:
        #        at=atoms_gen(item,at_ear,swap_dict)
        #        at_ear=at
            count_bonds_re=sro_count_bond(at,cutoff,pairs)
            delta,num_bond_disorder=func_all.order_para_cal(pairs,count_bonds_re,num_like,num_unlike,total_num_bonds) 
            print (pairs,count_bonds_re,delta)
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

        #with open('test.dat','a') as out:
        #    out.write(str (count_bonds))
        #    out.write('\t')
        #    out.write(str (atom_left))
        #    out.write('\n')
    else:
        print('cal_mode val not right')
        exit(1)
if __name__=="__main__":
    main()
