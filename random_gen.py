import numpy as np
from ase.neighborlist import neighbor_list
import random
import func_all
from ase import Atoms
from ase.io import write
import time

def tprofile_gen(Tmax,Tnum):
    '''

    '''
    '''
    mu_val = 2000
    sigma_val = 100
    mu, sigma = mu_val, sigma_val
    s = np.random.normal(mu_val, sigma_val, 1000)
    bins = np.linspace(0,2*mu_val,2*mu_val)
    gauss_vals = np.zeros(len(bins))
    gauss_vals = (1/(sigma * np.sqrt(2 * np.pi)) *\
            np.exp( - (bins - mu)**2 / (2 * sigma**2) ))
    max_gauss_val = np.max(gauss_vals)
    tfac = Tmax/max_gauss_val
    T_profile = np.zeros(len(bins))
    T_profile = tfac*gauss_vals
    max_tval = np.max(T_profile)
    for num,val in enumerate(T_profile):
        #an empirical condition to force the distribution to go to zero, 
        #if value is 2 order of magnitide lower than max amplitude of 
        #the gaussian distribution
        print (val,0.01*max_tval)
        if (val < (0.01*max_tval)):
            T_profile[num] = 0.
    '''
    alpha = 1
    T_profile = np.zeros(Tnum)
    for num in range(Tnum):
        T_profile[num] = Tmax/(1. + alpha*np.log(1+num))
    return T_profile
def swap_ideal(centre,elements,num_each,rand_pos,excluded_ele,atoms,elem_list,cutoff,pairs,bond_num,ideal_count,count_euc_ear,count_euc,T_profile,counter,anion_coord,anion):
    '''

    '''
    start1,end1 = func_all.start_end(elements,centre,num_each) 
    rand_swap = random.randint(0,len(excluded_ele)-1) #randomly choosing any of the excluded element, which will be swapped with other
    #finding the range for element, which will be swapped with 'other' element
    #print ('centre,rand_swap={},{}'.format(centre,rand_swap))
    start2,end2 = func_all.start_end(elements,excluded_ele[rand_swap],num_each)
    #print ('start1,end1={},{}'.format(start1,end1))
    at1 = random.randint(start1,end1-1)
    at2 = random.randint(start2,end2-1)
    #ensuring that at1 != at2 
    if (at1 == at2):
        while (at1 == at2):
            at2 = random.randint(start2,end2-1)
    else:
        pass 

    first = np.full(3,0,dtype=float)
    second = np.full(3,0,dtype=float)
    
    first = [a for a in rand_pos[at1]]
    second = [a for a in rand_pos[at2]]
    #print (first,second,pos[at1],pos[at2])
    for num,a in enumerate(second):
        rand_pos[at1][num] = a
    for num,a in enumerate(first):
        rand_pos[at2][num] = a
    #print ('at1,at2,pos[at1],pos[at2]={},{},{},{}'.format(at1,at2,pos[at1],pos[at2]))
    #print ('atoms.get_positions()[at1],atoms.get_positions()[at2]={},{}'.format(atoms.get_positions()[at1],atoms.get_positions()[at2]))
    at = Atoms(cell=atoms.get_cell()) 
    for a in range(len(rand_pos)):    
        at.extend(Atoms('{}'.format(elem_list[a]),positions=[(rand_pos[a][0],rand_pos[a][1],rand_pos[a][2])]))
    for a in range(len(anion_coord)):
        at.extend(Atoms('{}'.format(anion[0]),positions=[(anion_coord[a][0],anion_coord[a][1],anion_coord[a][2])]))
    at.set_pbc((True,True,True)) 
    #print ('at.get_positions()[at1],at.get_positions()[at2]={},{}'.format(at.get_positions()[at1],at.get_positions()[at2]))

    re_i,re_j=neighbor_list('ij',at,cutoff)
    #modify i and j to include only metals and remove anions
    re_i_mod = []
    re_j_mod = []
    for a,b in zip(re_i,re_j):
        if ( a < len(rand_pos) and b < len(rand_pos)):
            re_i_mod.append(a)
            re_j_mod.append(b)

    re_count_bond = func_all.bond_pair_count(pairs,re_i_mod,re_j_mod,elem_list)
    #print ('bond_num,re_count_bond,ideal_count={},{},{}'.format(bond_num,re_count_bond,ideal_count))
    
    #calculating difference between vectors
    dist_ear = np.linalg.norm(ideal_count-bond_num)
    dist_after = np.linalg.norm(ideal_count-re_count_bond)
    print ('dist_ear,dist_after,count_euc,count_euc_ear={},{},{},{}'.format(dist_ear,dist_after,count_euc,count_euc_ear))
    
    if (dist_ear >= dist_after or count_euc > (10*count_euc_ear)): #or count_euc > (10*count_euc_ear)
        return at,re_count_bond,dist_after,rand_pos
    elif (dist_ear < dist_after):
        #determining the probability
        #Temperature would be picked from the T_profile array, depending upon
        #the value of counter (or step number), if counter > len(T_profile), 
        #then T = 0
        if (counter >= len(T_profile)):
            #earlier swap is not accpeted and hence
            #positions are being changed to original.
            for num,a in enumerate(first):
                rand_pos[at1][num] = a
            for num,a in enumerate(second):
                rand_pos[at2][num] = a
            return atoms,bond_num,dist_ear,rand_pos
        else:
            temp = T_profile[counter]
            if (temp == 0):
                prob = 0. #to avoid zero division
            else:
                prob = np.exp((dist_ear - dist_after)/temp)
            accept = False
            max_num = 100
            prob_max = max_num * prob
            rand_val = random.randint(0,max_num)
            if (rand_val < int (prob_max)):
                accept = True
            else:
                accept = False
            if (accept == True):
                return at,re_count_bond,dist_after,rand_pos
            else:
                for num,a in enumerate(first):
                    rand_pos[at1][num] = a
                for num,a in enumerate(second):
                    rand_pos[at2][num] = a
                return atoms,bond_num,dist_ear,rand_pos

    else:
        #for num,a in enumerate(first):
        #    pos[at1][num] = a
        #for num,a in enumerate(second):
        #    pos[at2][num] = a
        #return atoms,bond_num,dist_ear,pos
        return at,re_count_bond,dist_after,rand_pos

def ideal(ideal_num,pairs,bond_num,nn_dict,num_each,elements,cn,rand_pos,atoms,elem_list,cutoff,Tmax,Tnum,anion_coord,anion):
    '''

    '''
    ideal_count = np.zeros(len(pairs))
    for num,item in enumerate(pairs):
        if (item.split('-')[0] == item.split('-')[1]):
            ideal_count[num] = ideal_num
        else:
            ideal_count[num] = 2*ideal_num
    #determination of array which contains the sequential decrease in T for simulated annealing
    
    T_profile = np.linspace(Tmax,0,Tnum)
    #T_profile = tprofile_gen(Tmax,Tnum)
    for item in T_profile:
        print (item)

    #T_profile = np.zeros(Tnum)
    euc_dist = np.linalg.norm(ideal_num - bond_num)
    start_euc = euc_dist
    #to track the steps spent at particular euler dist val
    euc_ear = euc_dist #intialisation
    count_euc_ear = 0 #steps spent at earlier euler dist val
    count_euc = 0 #steps spent at particular euler dist val
    counter = 0
    while (euc_dist != 0):
        #high_low = ideal_count - bond_num
        #print (high_low,ideal_count,bond_num)
        #inc = []
        #dec = []
        #inc_amt = []
        #dec_amt = []
        #for num,item in enumerate(high_low):
        #    if (item > 0):
        #        inc.append(pairs[num])
        #        inc_amt.append(item)
        #    elif (item < 0):
        #        dec.append(pairs[num])
        #        dec_amt.append(item)
        #print (inc,dec,inc_amt,dec_amt)
        #item_num = random.randint(0,len(dec)-1)
        #item = dec[item_num]
        #centre = item.split('-')[0]
        #other = item.split('-')[1]
        #print (centre,other)
        #finding bond instances which have 'other' element
        #and generating the excluded_ele = [] list
        #excluded_ele=[]
        #for val in inc:
        #    if (val.split('-')[0] == other):  #TODO other
        #        excluded_ele.append(val.split('-')[1])
        #    elif (val.split('-')[1] == other):
        #        excluded_ele.append(val.split('-')[0])
        #    else:
        #        pass
        #print ('exluded={}'.format(excluded_ele))
        #this condition is being added to deal with a situation,
        #when exluded_ele is empty
        #if (len(excluded_ele) == 0):
        #    el = random.randint(0,len(elements)-1)
        #    excluded_ele.append(elements[el])
        centre = elements[random.randint(0,len(elements)-1)]
        excluded_ele = []
        excluded_ele.append(elements[random.randint(0,len(elements)-1)])
        at,re_bond_num,euc_dist,pos = swap_ideal(centre,elements,num_each,rand_pos,\
                excluded_ele,atoms,elem_list,cutoff,pairs,bond_num,ideal_count,count_euc_ear,count_euc,T_profile,counter,anion_coord,anion)
        counter += 1
        #if euclidian distance is unchanged
        if (euc_dist == euc_ear):
            count_euc += 1
        #if euclidian distance changes
        if (euc_dist != euc_ear):
            count_euc_ear = count_euc
            count_euc = 0
        print ('euc_dist,euc_ear,count_euc,count_euc_ear={},{},{},{}'.format(euc_dist,euc_ear,count_euc,count_euc_ear))
        euc_ear = euc_dist
        #changing
        atoms = at
        bond_num = re_bond_num
        rand_pos = pos
        with open ('euc_distance.csv','a') as inp:
            inp.write(str(euc_dist))
            inp.write('\n')
    return atoms,bond_num,start_euc,counter
        #ele_id = func_all.ele_indx(elements,centre)
        #start = 0
        #end = 0
        #for a in range(len(elements)):
        #    if (a < ele_id):
        #        start = start + num_each[a]
        #        end = end + num_each[a]
        #    elif (a == ele_id):
        #        end = end + num_each[a]
        #    else:
        #        pass
        #start,end = func_all.start_end(elements,centre,num_each)
        #print ('centre,ele_id,num_each,start,end={},{},{},{},{}'.format(centre,ele_id,num_each,start,end))
        #centre_val = []
        #other_val = []
        #for a in range(start,end):
        #    nnlist1=np.array(range(cn),dtype='int')
        #    for b,c in enumerate(nn_dict[a]):
        #        nnlist1[b] = c 
        #    for d in nnlist1:
        #        if (list_element[d] == other):
        #            rand_swap = random.randint(0,len(excluded_ele)) #randomly choosing any of the excluded element, which will be swapped with other
                    #finding the range for element, which will be swapped with 'other' element
        #            start2,end2 = func_all.start_end(elements,exluded_ele[rand_swap],num_each)
        #            other_swap = random.randint(start2,end2) #element entry which will be swapped with 'other' element or 'd'
                    #ensuring that other_swap != a 
        #            if ( other_swap == a):
        #                while (other_swap == a):
        #                    other_swap = random.randint(start2,end2)
        #            else:
        #                pass 
                    
                    #nnlist2=np.array(range(cn),dtype='int')
                    #for m,n in enumerate(nn_dict[d]):
                    #    nnlist2[m] = n
                    #checking, if any ele in nnlist2 is excluded_ele
                    #count = 0
                    #for m in nnlist2:
                    #    for n in excluded_ele:
                    #        if (n == list_element[m]):
                    #            break
                    #        else:
                    #            count += 1
                    #print ('count={}'.format(count))
                    #if ( count == ((cn * len(excluded_ele)) - 1)):         
                    #    centre_val.append(a)
                    #    other_val.append(d)
        #print ('central_val,other_val={},{}'.format(centre_val,other_val))
            
def roulette(atom_left,elements):
    '''

    '''
    max_num = 10000
    #fraction of atom_left
    frac=np.zeros(len(elements))
    tot=sum(a for a in atom_left)
    for num,item in enumerate(atom_left):
        frac[num] = item/tot
    add_frac=np.zeros(len(elements))
    count=0.
    for num,item in enumerate(frac):
        add_frac[num] = frac[num] + count
        count=add_frac[num]
    val_range=np.zeros(len(elements))
    for num,item in enumerate(add_frac):
        val_range[num] = int  (max_num * item)


    rand_val = random.randint(0,max_num)
    at = -1
    for num,item in enumerate(val_range):
        if (item >= rand_val):
            at = num
            break
        else:
            pass
    #print (frac,add_frac,val_range,rand_val,at)
    return at

def check_bond_count(list_element,i,j,pairs):
    temp_ele=set(list_element)
    #since set object is not subscriptable, we need to generate a list
    uniq_ele=[]
    for g in temp_ele:
        uniq_ele.append(g)
    #determination of pairs with elements added so far and proxy element
    check_num_pairs=func_all.num_det(uniq_ele)
    check_num_like=len(uniq_ele)
    check_num_unlike=check_num_pairs-check_num_like
    check_pairs=[]
    #generating bonding pairs
    check_pairs=func_all.pairs_gen(uniq_ele,check_num_pairs)
    check_count_bonds=[]
    #start8=time.time() #TODO remove
    check_count_bond = func_all.bond_pair_count(check_pairs,i,j,\
            list_element) #check_i -> i and check_j -> j
    #end8=time.time()
    #print ('time-taken8={}'.format(end8-start8))
    check_indx=[]
    check_bond_num=[]
    for g,h in enumerate(check_pairs):
        temp1=h.split('-')[0]
        temp2=h.split('-')[1]
        for k,m in enumerate(pairs):
            temp3=m.split('-')[0]
            temp4=m.split('-')[1]
            if (temp3==temp1 and temp4==temp2) or (temp3==temp2 and temp4==temp1):
                check_indx.append(k)
                check_bond_num.append(check_count_bond[g])
                break
            else:
                pass
    #checking if the bond_number is within the max bond number limit
    #max bond limit check should only be done for bond formed due to elem_id and elem_id2
    #pindex=pair_index(pairs,elements[elem_id],elements[elem_id2]) #det pair index in pairs[]
    #print('check_indx,check_bond_num,ideal_num_counter={},{},{}'.format(\
    #        check_indx,check_bond_num,ideal_num_counter))
    return check_indx,check_bond_num
def check_bond_avl(list_element,i,j,elements,pos_change,elem_id,ideal_num_counter,pairs,atom_left):
    ''' 
    Input
    pos_change: position in list_element, where proxy_ele
    has been replaced with the element[elem_id]
    Returns:
    bond_avl: Boolean parameter. True: available and 
                                 False: not available 
    '''
    #if other elements are not left in the atom_left array
    ele_unavail=0
    for num,item in enumerate(atom_left):
        if(num == elem_id):
            continue        #because it has been ensured in random_gen function that elem_id is available.
        else:
            if(item == 0):
                ele_unavail+=1
            else:
                pass
    if(ele_unavail > len(elements)-2):    #no other elements are left 
        bond_avl=True
        return bond_avl
    else: #this condition is only executed, when other elements are available.
        #start7=time.time() #TODO remove
        check_indx,check_bond_num=check_bond_count(list_element,i,j,pairs)
        #end7=time.time()
        #print ('time-taken7={}'.format(end7-start7))
        for k,m in enumerate(check_indx):
            if(check_bond_num[k] > ideal_num_counter[m]):
                #print('In check_bond_avl fun:k,check_bond_num[k],m,ideal_num_counter[m],elem_id,elements={},{},{},{},{},{}'.format(k\
                #        ,check_bond_num[k],m,ideal_num_counter[m],elem_id,elements))
                #checking, whether intro of other elements at the lattice also leads to the 
                #increase in the bond-count beyond max limit, in that case, each element, except
                #the earlier element would be introduced and bond-count would be checked,
                #If, all other bond-counts are also going to increase beyond max-value,
                #then, earlier element would be chosen.
                element_counter=0 #counter to store the instances, when bond-val is > max-limit
                for a in range(len(elements)):
                    if(a==elem_id):
                        pass
                    else:
                        list_element[pos_change]=elements[a]
                        check_indx_new,check_bond_new=check_bond_count(list_element,i,j,pairs)
                        #print('element[a],check_bond_new,check_indx_new={},{},{}'.format(elements[a],\
                        #        check_bond_new,check_indx_new))
                        for x,y in enumerate(check_indx_new):
                            #If the bond-val is > max_limit or there is no atom left (to avoid recussive loop)
                            if(check_bond_new[x] > ideal_num_counter[y] or atom_left[a]==0): 
                                element_counter+=1
                                break
                #print('element_counter={}'.format(element_counter))
                list_element[pos_change]=elements[elem_id] #undoing the change in list_element[]
                if(element_counter==len(elements)-1):
                    #print('element_counter reached max value={}'.format(element_counter))
                    bond_avl=True
                    #list_element[pos_change]=elements[elem_id] #TODO check, if this needs to be uncommented!
                else:
                    bond_avl=False
                break 
            else:
                bond_avl=True
        return bond_avl
def arrange_atoms(elements,num_each,list_element,check_pos,atoms):
    re_list_element=np.array(range(len(atoms)),dtype='str')
    atom_re=Atoms(cell=atoms.get_cell())
    counter=0
    for a in range(len(elements)):
        for b in range(len(atoms)):
            if(list_element[b] == elements[a]):
                atom_re.extend(Atoms('{}'.format(list_element[b]),\
                        positions=[(check_pos[b][0],check_pos[b][1],\
                        check_pos[b][2])]))
                re_list_element[counter]=list_element[b]
                counter+=1
            else:
                pass
    atom_re.set_pbc(111)
    return atom_re,re_list_element    
def random_struc_gen(proxy_ele,elements,num_each,at,cutoff,pairs,num_unlike,num_like,cn,Tmax,Tnum,i,j,elem_list,rand_pos,anion_coord,anion):
     
    #option-1: simply shuffling the atomic positions to generate the initial random structure
    #determining number of each chemical species
    #num_each = func_all.number_det(at,elements)

    #list of atoms in a sequence
    #elem_list = func_all.element_assign(at,elements,num_each)

    #rand_pos = func_all.random_gen(at)

    #elem_id = func_all.elem_dict(elements)

    #atoms = Atoms(cell=at.get_cell())
    #for a in range(len(rand_pos)):
    #    atoms.extend(Atoms('{}'.format(elem_list[a]),positions=[(rand_pos[a][0],rand_pos[a][1],rand_pos[a][2])]))
    #atoms.set_pbc((True,True,True))
    #print (atoms)
    #atoms.set_pbc(111)
    #i,j=neighbor_list('ij',atoms,cutoff) 
    ideal_num=(len(i))/((2*num_unlike)+num_like) #ideal no of like bonds for complete disorder
    #checking, if the ideal_num is an integer, if not program need to be terminated and
    #different inp.cfg file need to be generated, such that total number of bonds must 
    #be divisible with 2*k1+k2
    check_int = len(i) % ((2*num_unlike)+num_like)
    if (check_int != 0):
        print ('Make sure that inp.cfg contains no. of atoms such that number of bonds should be divisible with {}'.format(2*num_unlike+num_like))
        exit(1)
    count_bond = func_all.bond_pair_count(pairs,i,j,elem_list)
    nn_dict=func_all.nn_dict_det(i,j)
    atoms_final,count_bond_final,euc_start,counter_reach= ideal(ideal_num,pairs,count_bond,nn_dict,num_each,elements,cn,rand_pos,at,elem_list,cutoff,Tmax,Tnum,anion_coord,anion)
    #with open('num_steps_req.dat','a') as inp:
    #    inp.write(str(euc_start))
    #    inp.write('\t')
    #    inp.write(str(counter_reach))
    #    inp.write('\n')

    atom_left = np.zeros(len(num_each)) #no meaning, it is here simply to avoid error
    delta,num_bond_disorder=func_all.order_para_cal(pairs,count_bond_final,\
            num_like,num_unlike,len(i))
    if(delta == 0.0):
        write('{}_random.xyz'.format(delta),atoms_final)
        write('{}_random.vasp'.format(delta),atoms_final)

    '''
    #option-2: engineered shuffling
    
    #TEST: This is being included to ensure that proxy element might not be in the original list
    for item in elements:
        if(item == proxy_ele):
            print("Change the proxy element")
            exit(1)
        else:
            pass
    #start1=time.time() #TODO remove this
    positions=at.get_positions()
    list_element=np.array(range(len(at)),dtype='str')
    for num in range(len(at)):
        list_element[num]=proxy_ele
    #print(list_element)
    #1/0
    atom_left = np.zeros(len(num_each))
    for num,item in enumerate(num_each):
        atom_left[num] = item
    i,j=neighbor_list('ij',at,cutoff)
    #end1=time.time() #TODO remove this
    #print('time taken1={}',format(end1-start1))
    ideal_num=(len(i))/((2*num_unlike)+num_like) #ideal no of like bonds for complete disorder
    print ('ideal_num,len(i),num_like,num_unlike = {},{},{},{}'.format(ideal_num,len(i),num_like,num_unlike))
    #ideal bond counter for all bonds (in complete disorder)
    ideal_num_counter=np.zeros(num_like+num_unlike)
    for num,item in enumerate(pairs):
        temp1=item.split('-')[0]
        temp2=item.split('-')[1]
        if(temp1 == temp2):
            ideal_num_counter[num]=ideal_num
        else:
            ideal_num_counter[num]=2*ideal_num
    #initialising the bond_counter, which will the store the bond_instances
    #bond_counter=np.zeros(len(pairs))
    nn_dict=func_all.nn_dict_det(i,j)
    forbidden_list=np.full(len(at),-1,dtype=int)
    counter_fl=0
    for a in range(len(at)):
        #checking if 'a' is already in the forbidden_list (check_forbd=True and False otherwise)
        check_forbd=func_all.fl_check(forbidden_list,a)
        #print(a,check_forbd,func_all.forbidden_list_count(forbidden_list),atom_left)
        if(check_forbd == True):
            continue
        elif(check_forbd == False):
            #print ('a,check_forbd,fl={},{},{}'.format(a,check_forbd,forbidden_list))
            #elem_id=random.randint(0,len(elements)-1)
            #using dynamic roulette wheel for selecting element below
            elem_id = roulette(atom_left,elements)
            #putting a condition that if one of the element is zero, 
            #that particular atom needs to be avoided.
            if(atom_left[elem_id]==0):
                while (atom_left[elem_id] == 0):
                    #elem_id=random.randint(0,len(elements)-1) #ensuring that atom chosen is available
                    #using dynamic roulette wheel for selecting element below
                    elem_id = roulette(atom_left,elements)
                    if(atom_left[elem_id] > 0):
                        break
                    else:
                        pass
            else:
                pass
            #print('atom_left={}'.format(atom_left)) #TODO remove
            #start2=time.time() #TODO remove this
            for b,c in enumerate(atom_left):
                if(b == elem_id and c > 0):
                    #print('b,c={},{}'.format(b,c)) #TODO comment
                    #it is being assumed that number of nn is CN
                    nnlist=np.array(range(cn),dtype='int')
                    for num,item in enumerate(nn_dict[a]):
                        nnlist[num]=item
                    #nnlist=nn_dict[a] This is being removed due to anomaly in removing elements
                    #print('before putting central atom={}'.format(list_element[a]))
                    list_element[a]=elements[elem_id] #putting required element at the particular place
                    #print('after putting central atom={}'.format(list_element[a]))
                    #increasing the counter for the central atom
                    atom_left[b] = atom_left[b]-1
                    #print('nnlist={}'.format(nnlist))
                    #adding the central atom to the forbidden_list
                    forbidden_list[counter_fl]=a
                    counter_fl+=1
                    to_remove=[] #list to store the indices of the entries to be removed from nnlist
                    #start3=time.time() #TODO remove 
                    for num,item in enumerate(nnlist):#if items in nnlist are already in forbidden_list
                        check_forbd2=func_all.fl_check(forbidden_list,item)
                        if(check_forbd2 == True):
                            #print('nnlist item to be removed={}'.format(item))
                            #nnlist.pop(num)
                            to_remove.append(num)
                            continue
                        elif(check_forbd2 == False):
                            pass
                        else:
                            print('check_forbd2 has an unexpected value')
                            exit(1)
                    nnlist_rem=nnlist[~np.isin(np.arange(nnlist.size), to_remove)]
                    #end3=time.time() #TODO remove
                    #print('time-taken3={}'.format(end3-start3))
                    #print('len of mod nnlist={},{}'.format(len(nnlist_rem),nnlist_rem))
                    
                    #deleting the duplicate entries in nnlist_mod 
                    nnlist_mod = np.unique(nnlist_rem)
                    #print('len of mod nnlist={},{}'.format(len(nnlist_mod),nnlist_mod))
                    #TODO add the condition that if nnlist length is zero and due to introducing
                    #particular atom, bond-instances increases beyond the specified limit, then 
                    #central atom should be reverted back to the proxy element and atom_left counter
                    #need to be increased.
                    if(len(nnlist_mod) == 0):
                        bond_avl=check_bond_avl(list_element,i,j,elements,a,elem_id,\
                                ideal_num_counter,pairs,atom_left)
                        #print ('bond_avl before={}'.format(bond_avl)) #TODO remove this
                        if(bond_avl==False):
                            list_element[a]=proxy_ele
                            forbidden_list[counter_fl-1]=-1
                            counter_fl-=1
                            atom_left[b] = atom_left[b]+1
                        else:
                            pass
                        #print('list_element[a],atom_left[b],counter_fl,fl[cfl]={},{},{},{}'.format(list_element[a],\
                        #        atom_left[b],counter_fl,forbidden_list[counter_fl])) 
                    else:
                        pass
                    #start4=time.time() #TODO remove 
                    for num,item in enumerate(nnlist_mod):
                        #print(num,item)
                        #checking if the element to be put in NN
                        #is available and bond_counter allows the bond
                        atom_avl=False
                        bond_avl=False
                        while (atom_avl == False and len(nnlist_mod) > 0) or (bond_avl == False and len(nnlist_mod) > 0):
                            #start5=time.time() #TODO remove
                            #TODO remove this after the check before elem_id2 line
                            #temp_ele=set(list_element)
                            #since set object is not subscriptable, we need to generate a list
                            #uniq_ele=[]
                            #for g in temp_ele:
                            #    uniq_ele.append(g)
                            #determination of pairs with elements added so far and proxy element
                            #check_num_pairs=func_all.num_det(uniq_ele)
                            #check_num_like=len(uniq_ele)
                            #check_num_unlike=check_num_pairs-check_num_like
                            #check_pairs=[]
                            #generating bonding pairs
                            #check_pairs=func_all.pairs_gen(uniq_ele,check_num_pairs)
                            #check_count_bonds=[]
                            #check_count_bond = func_all.bond_pair_count(check_pairs,i,j,\
                            #        list_element) #check_i -> i and check_j -> j
                            #print('before introducing the atom in NN={},{}'.format(check_pairs,check_count_bond))
                            #TODO remove above
                            
                            #elem_id2=random.randint(0,len(elements)-1)
                            #using dynamics roulette wheel for selecting element
                            elem_id2 = roulette(atom_left,elements)
                            
                            #el_indx=indx_find(num_each,item,elements) #finding indx for atom_left
                            if(atom_left[elem_id2] > 0):
                                atom_avl=True
                                #TODO remove this after check
                                #proxy_count=func_all.count_proxy(proxy_ele,list_element)
                                #print('before swap, proxy count & list_element[item]={},{}'.format(proxy_count,\
                                #        list_element[item]))
                                list_element[item]=elements[elem_id2] #elem is being put at NN site as trial.
                                #TODO remove this after check
                                #proxy_count=func_all.count_proxy(proxy_ele,list_element)
                                #print('after swap, proxy count & list_element[item]={},{}'.format(proxy_count,\
                                #        list_element[item]))
                                #print('elem_id2,atom_left[elem_id2],item,list_element[item]={},{},{},{}'.format\
                                #        (elem_id2,atom_left[elem_id2],item,list_element[item]))
                            else:
                                atom_avl=False
                                #print('No atom_left={}'.format(atom_left[elem_id2]))
                                continue #other element need to be chosen
                            #start6=time.time() #TODO remove
                            bond_avl=check_bond_avl(list_element,i,j,elements,\
                                    item,elem_id2,ideal_num_counter,pairs,atom_left)
                            #print ('bond_avl after={}'.format(bond_avl)) #TODO remove this
                            #end6=time.time()
                            #print('time-taken6={}'.format(end6-start6))
                            #print('atom_avl,bond_avl={},{}'.format(atom_avl,bond_avl))
                            #print('list_element[item] after bond_avl check={}'.format(list_element[item]))
                            #If, adding an atom to NN site increases the bond count beyond the max no.
                            if(bond_avl == False):
                                #print('before reverting={}'.format(list_element[item]))
                                list_element[item]=proxy_ele #proxy ele is again put in the place of trail atom
                                #print('after reverting={}'.format(list_element[item]))
                            elif(bond_avl == True):
                                atom_left[elem_id2]=atom_left[elem_id2]-1
                                forbidden_list[counter_fl]=item
                                counter_fl+=1
                            #end5=time.time()
                            #print('time-taken5={}'.format(end5-start5))
                    #end4=time.time()
                    #print('time-taken4={}'.format(end4-start4))
                                #print('forbidden_list len after swap,atom_left={},{}'.format(\
                                #        func_all.forbidden_list_count(forbidden_list),atom_left))
                            #adding the NN element to the forbidden_list
    #TO deal with a bug in the system, where some proxy element are left, even going  
    #the earlier loop, proxy elements would be removed and replaced with the atom_left
    #atom_left[0]=5
    #list_element[50]='Ca'
    #list_element[67]='Ca'
    #list_element[3]='Ca'
    #list_element[45]='Ca'
    #list_element[78]='Ca'
            #end2=time.time() #TODO remove this
            #print('time-taken2={}'.format(end2-start2))
    atom_left_counter=sum(item for item in atom_left)
    if(atom_left_counter > 0):
        #find the places where proxy_ele are present
        proxy_place=[]
        proxy_place=[num for num,item in enumerate(list_element) if item==proxy_ele]
        #print('proxy_place,len(proxy_place),atom_left_counter={},{},{}'.format(proxy_place,\
        #        len(proxy_place),atom_left_counter))
        if(len(proxy_place) != atom_left_counter):
            print('Issue with proxy_place')
            exit(1)
        else:
            pass
        left_counter=0
        while (left_counter < atom_left_counter):
            el_id=random.randint(0,len(elements)-1)
            if(atom_left[el_id] > 0):
                list_element[proxy_place[left_counter]]=elements[el_id]
                left_counter+=1
    else:
        pass
    #TEST
    temp_ele=set(list_element)
    #since set object is not subscriptable, we need to generate a list
    uniq_ele=[]
    for g in temp_ele:
        uniq_ele.append(g)
    check_at = Atoms(cell=at.get_cell())
    check_pos = []
    for g in range(len(positions)):
        check_pos.append(positions[g])

    for g in range(len(check_pos)):
        check_at.extend(Atoms('{}'.format(list_element[g]),\
                positions=[(check_pos[g][0],check_pos[g][1],\
                check_pos[g][2])]))
    check_at.set_pbc(111)
    check_i,check_j=neighbor_list('ij',check_at,cutoff)
    
    #determination of pairs with elements added so far and proxy element
    check_num_pairs=func_all.num_det(uniq_ele)
    check_num_like=len(uniq_ele)
    check_num_unlike=check_num_pairs-check_num_like
    check_pairs=[]
    #generating bonding pairs
    check_pairs=func_all.pairs_gen(uniq_ele,check_num_pairs)
    check_count_bonds=[]
    check_count_bond = func_all.bond_pair_count(check_pairs,check_i,check_j,list_element)
    atom_re,re_list_element=arrange_atoms(elements,num_each,list_element,check_pos,check_at)
    
    #check
    re_i,re_j=neighbor_list('ij',atom_re,cutoff)
    re_count_bond=func_all.bond_pair_count(check_pairs,re_i,re_j,re_list_element)
    #print(check_pairs,check_count_bond,atom_left)
    #print('modified list_element bond count={}'.format(re_count_bond))
    
    atoms_final,count_bond_final,counter_reach = ideal(ideal_num,pairs,check_count_bond,nn_dict,num_each,elements,cn,check_pos,check_at,list_element,cutoff)
    #writting the number of steps to reach delta = 0
    with open('num_steps_req.dat','a') as inp:
        inp.write(str(counter_reach))
        inp.write('\n')
    #writing atoms object with its order parameter in the file
    delta,num_bond_disorder=func_all.order_para_cal(pairs,count_bond_final,\
            num_like,num_unlike,len(re_i))
    #note that num_bond_disorder is n parameter in delta parameter.
    #print (delta)
    if(delta == 0.0):
        write('{}_random.xyz'.format(delta),atom_re)
        write('{}_random.vasp'.format(delta),atom_re)
    #print ('returning val = {},{},{},{},{}'.format(check_count_bond,atom_left,delta,ideal_num,num_each))
    #ideal(ideal_num,pairs,check_count_bond,nn_dict,num_each,elements,cn,check_pos,check_at,list_element,cutoff)
    '''
    return count_bond_final,atom_left,delta,ideal_num
    
