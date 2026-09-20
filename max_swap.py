import numpy as np
import random
import func_all

def CA_SP_det(bond1,bond2,pref_pair):
    central_pair = []
    swap_pair = []
    first_x = bond1.split('-')[0]
    second_x = bond1.split('-')[1]
    if(first_x == pref_pair[0] or first_x == pref_pair[1]):
        unpref_at1 = second_x
        pref_at1 = first_x
    elif(second_x == pref_pair[0] or second_x == pref_pair[1]):
        pref_at1 = second_x
        unpref_at1 = first_x
    else:
        raise Exception('Issue with the unpref_at1 determination')

    first_y = bond2.split('-')[0]
    second_y = bond2.split('-')[1]
    if(first_y == pref_pair[0] or first_y == pref_pair[1]):
        pref_at2=first_y
        unpref_at2=second_y
    elif(second_y == pref_pair[0] or second_y == pref_pair[1]):
        pref_at2=second_y
        unpref_at2=first_y
    else:
        raise Exception('Issue with pref_at2 determination')
    if(pref_pair[0] != pref_pair[1]):
        swap_pair.append(unpref_at1)
        swap_pair.append(pref_at2)
        central_pair.append(pref_at1)
        central_pair.append(unpref_at2)
    elif(pref_pair[0] == pref_pair[1]):
        swap_pair.append(unpref_at1)
        swap_pair.append(pref_at1)
        central_pair.append(pref_at1)
        central_pair.append(unpref_at2)
    else:
        raise Exception('Issue with CA and SP determination')

    return central_pair,swap_pair
def bond_swap_metric(group1,group2,pref_pair):
    '''
    Parameters:
        group1:
        group2:
        pref_pair_req:
    Returns:
        group1_new:
        group2_new:
    '''
    if (pref_pair[0] == pref_pair[1]):
        like=True
        unlike=False
    elif(pref_pair[0] != pref_pair[1]):
        like=False
        unlike=True
    else:
        raise Exception('Issue with pref_pair list')

    metric=np.empty(len(group1)*len(group2),dtype='int')
    spairs=[]
    group1_new = []
    group2_new = []
    count = 0
    for num1,item1 in enumerate(group1):
        for num2,item2 in enumerate(group2):
            central_pair,swap_pair = CA_SP_det(item1,item2,pref_pair)
            spair = []
            '''
            central_pair[0] - swap_pair[1] forms the preferred bonds and there is associated rise in the propensity
            of central_pair[1] - swap_pait[0] bond, with decrease in the central_pair[0] - swap_pair [0] and 
            central_pair[1] - swap_pair[1] bonds.
            '''
            if (like == True and unlike == False):
                '''
                for calculation of metric, we assign -1 to the undesired bond change (i.e., increase in the unlike
                bond), while +1 to the desired change (increase in the like bond) and +2 for the increase in the 
                preferred bond number
                '''
                if (central_pair[0] == swap_pair[0]):
                    first = -1 #since this desired type of the bond is going to be decreased
                if (central_pair[0] != swap_pair[0]):
                    first = 1 #since this undesired type of the bond is going to be decreased
                if (central_pair[1] == swap_pair[1]):
                    second = -1 #since this desired type of the bond is going to be decreased
                if (central_pair[1] != swap_pair[1]):
                    second = 1 #since this undesired type of the bond is going to be decreased
                if (central_pair[0] == swap_pair[1]):
                    third = 1 + 1#since this desired type of the bond (&preferred bond) is going to be increased
                if (central_pair[0] != swap_pair[1]):
                    third = -1 #since this desired type of the bond is going to be increased, but hopefully it wouldn't happen!
                if (central_pair[1] == swap_pair[0]):
                    fourth = 1 #since this desired type of the bond is going to be increased.
                if (central_pair[1] != swap_pair[0]):
                    fourth = -1 #since this undesired type of the bond is going to be increased.
            
            if (like == False and unlike == True):
                '''
                for calculation of metric, we assign -1 to the undesired bond change (i.e., increase in the like
                bond), while +1 to the desired change (increase in the unlike bond) and +2 for the increase in the 
                preferred bond number.
                '''
                if (central_pair[0] == swap_pair[0]):
                    first = 1 #since this undesired type of the bond is going to be decreased
                if (central_pair[0] != swap_pair[0]):
                    first = -1 #since this desired type of the bond is going to be decreased
                if (central_pair[1] == swap_pair[1]):
                    second = 1 #since this undesired type of the bond is going to be decreased
                if (central_pair[1] != swap_pair[1]):
                    second = -1 #since this desired type of the bond is going to be decreased
                if (central_pair[0] == swap_pair[1]):
                    third = -1 #since this desired type of the bond is going to be increased, but hopefully it wouldn't happen!
                if (central_pair[0] != swap_pair[1]):
                    third = 1 + 1 #since this desired type of the bond (&preferred bond) is going to be increased.
                if (central_pair[1] == swap_pair[0]):
                    fourth = -1 #since this undesired type of the bond is going to be increased.
                if (central_pair[1] != swap_pair[0]):
                    if ((central_pair[1] == pref_pair[0] and swap_pair[0] == pref_pair[1]) or (central_pair[1] == pref_pair[1] and swap_pair[0] == pref_pair[0])): 
                        fourth = 1 + 1  #since this desired type of the bond (& preferred bond)  is going to be increased.
                    else:
                        fourth = 1
            
            spair.append(num1)
            spair.append(num2)
            metric[count] = first + second + third + fourth
            count += 1
            
            spairs.append(spair)
    

    #sorting the metric and arranging bonds 
    sorted_index_pos = np.flip([index for index, num in sorted(enumerate(metric), key=lambda x: x[-1])])
    #print (sorted_index_pos)
    for val in sorted_index_pos:
        group1_new.append(group1[spairs[val][0]])
        group2_new.append(group2[spairs[val][1]])
    return group1_new,group2_new 

def group_pair_metric(bond1, bond2, pref_pair):
    central_pair, swap_pair = CA_SP_det(bond1, bond2, pref_pair)
    if pref_pair[0] == pref_pair[1]:
        if central_pair[0] == swap_pair[0]:
            first = -1
        else:
            first = 1
        if central_pair[1] == swap_pair[1]:
            second = -1
        else:
            second = 1
        if central_pair[0] == swap_pair[1]:
            third = 2
        else:
            third = -1
        if central_pair[1] == swap_pair[0]:
            fourth = 1
        else:
            fourth = -1
    elif pref_pair[0] != pref_pair[1]:
        if central_pair[0] == swap_pair[0]:
            first = 1
        else:
            first = -1
        if central_pair[1] == swap_pair[1]:
            second = 1
        else:
            second = -1
        if central_pair[0] == swap_pair[1]:
            third = -1
        else:
            third = 2
        if central_pair[1] == swap_pair[0]:
            fourth = -1
        else:
            if ((central_pair[1] == pref_pair[0] and swap_pair[0] == pref_pair[1]) or
                (central_pair[1] == pref_pair[1] and swap_pair[0] == pref_pair[0])):
                fourth = 2
            else:
                fourth = 1
    else:
        raise Exception('Issue with pref_pair list')
    metric = first + second + third + fourth
    return metric
def group_determination(pref_pair_req,pairs):
    '''
    Parameters:
        pref_pair_req:
        pairs:
    
    Returns: 
        group1:
        group2:
        pref_pair:
    
    '''
    pref_pair=[]
    #ensuring that pref_pair has same order as in pair[]
    for val in pairs:
        val_temp1=val.split('-')[0]
        val_temp2=val.split('-')[1]
        if pref_pair_req[0] == val_temp1 and pref_pair_req[1] == val_temp2:
            pref_pair.append(pref_pair_req[0])
            pref_pair.append(pref_pair_req[1])
        elif pref_pair_req[0] == val_temp2 and pref_pair_req[1] == val_temp1:
            pref_pair.append(pref_pair_req[1])
            pref_pair.append(pref_pair_req[0])
        else:
            pass
    #TEST
    if(len(pref_pair) == 0):
        raise Exception('Issue with pref_pair list') 
    
    group1_tmp=[]
    group2_tmp=[]
    count_like=0
    #group1 and group2 are being generated. group1 and group2 represent the bonds, between which 
    #swap would be carried out to increase the preferred bonds
    for item in pairs:
        if(pref_pair[0] != pref_pair[1]):
            temp1=item.split('-')[0]
            temp2=item.split('-')[1]
            if(temp1 == pref_pair[0] or temp2 == pref_pair[0]): #pref_pair[0] --> group1
                if((temp1 == pref_pair[0] and temp2 == pref_pair[1]) or (temp1==pref_pair[1] and temp2==pref_pair[0])):
                    pass #pref_pair is not added in the group1
                else:
                    group1_tmp.append(item)
            if(temp1 == pref_pair[1] or temp2 == pref_pair[1]): #pref_pair[1] --> group2
                if((temp1 == pref_pair[0] and temp2 == pref_pair[1]) or (temp1==pref_pair[1] and temp2==pref_pair[0])):
                    pass #pref_pair is not added in the group1
                else:
                    group2_tmp.append(item)
        if(pref_pair[0] == pref_pair[1]):
            temp1=item.split('-')[0]
            temp2=item.split('-')[1]
            if(temp1==pref_pair[0] or temp2==pref_pair[1]):
                if(temp1==pref_pair[0] and temp2==pref_pair[1]):
                    pass
                elif(count_like%2 != 0):
                    group1_tmp.append(item)
                    count_like+=1
                elif(count_like%2 == 0):
                    group2_tmp.append(item)
                    count_like+=1
    
    group1,group2 = bond_swap_metric(group1_tmp,group2_tmp,pref_pair)
    #print ('group1, group2={}{}'.format(group1,group2))    
    #TODO possibly the below-given code is redundant as we are simply
    #generating the group1 and group2 with an idea that swap between
    #bonds are carried out to increase the value of a matric which
    #is decided on the basis of sign of delta parameter desired. Also
    #it is seen that generally for constant compositon cases, swap 
    #between one pair if bond is good enough to fill the forbidden list!

    #It is being ensured here that the length of group1 and group2
    #are equal, it is particularly important, when the preferred
    #bond is a 'like' bond. As the odd number of bond-pairs are 
    #possible for generation of like bonds, if the number of 
    #elements in the system are even 
    if(len(group1) != len(group2)):
        temp=-1
        num=-1
        if(len(group1) > len(group2)):
            temp=len(group1)-len(group2)
            #temp should always be equal to 1
            #as to make any number even from odd only one should be removed
            if (temp != 1):
                print ('Issue with group1 and group2 generation')
                exit(1)
            else:
                pass
            num=random.randint(0,(len(group1)-1))
            group1.pop(num)
             
        elif(len(group1) < len(group2)):
            temp=len(group2)-len(group1)
            #temp should always be equal to 1
            #as to make any number even from odd only one should be removed
            if (temp != 1):
                raise Exception('Issue with group1 and group2 generation-2')
            else:
                pass
            num=random.randint(0,(len(group2)-1))
            group2.pop(num)

    return group1,group2,pref_pair

# group_determination() and group_pair_is_useful() stay exactly as you have them —
# not touched.
def group_pair_is_useful(bond1, bond2, pref_pair):
    """
    Determine whether a group1-group2 bond combination is
    worth sending to the expensive atom-level search.

    Returns:
        True  -> perform atom-level search
        False -> reject immediately
    """

    metric = group_pair_metric(
        bond1,
        bond2,
        pref_pair
    )

    if metric > 0:
        return True
    else:
        return False
 



# CA_SP_det, bond_swap_metric, group_pair_metric, group_determination,
# group_pair_is_useful all stay exactly as you have them - not touched.

def swap_pos(pref_pair_req,pairs,elements,num_each,nn_dict,count_bonds):
    total_CA1_selected = 0
    successful_swaps = 0
    failed_swaps = 0
    terminated_due_to_full = 0
    rows=len(elements)
    cols=[]
    for item in num_each:
        cols.append(item)
    forbidden_vals=[[-1 for a in range(cols[b])] for b in range(rows)]
    counts=np.zeros(len(elements))
    forbidden_list=[]
    num_swap=0
    num_bond_mod=np.zeros(len(pairs))
    group1,group2,pref_pair = group_determination(pref_pair_req,pairs)
    swap_dict={}
    DEBUG = False

    total_pairs_considered = 0
    pairs_with_swap = 0
    pairs_without_swap = 0

    indx_red1 = -1
    indx_red2 = -1
    indx_inc1 = -1
    indx_inc2 = -1

    # ---------------------------------------------------------------
    # group1[k] / group2[k] are a MATCHED pair - group_determination()
    # (via bond_swap_metric) already built and sorted them TOGETHER,
    # best metric first. Walk that list in order. Do NOT cross-combine
    # group1[num_x] with group2[num_y] for num_x != num_y - that
    # breaks the pairing the sort was built around and is what was
    # causing the wrong bonds to move.
    # ---------------------------------------------------------------
    for k in range(len(group1)):

        g1 = group1[k]
        g2 = group2[k]

        total_pairs_considered += 1

        pair_metric = group_pair_metric(g1, g2, pref_pair)

        if pair_metric <= 0:
            # sorted descending -> everything from here on is also <= 0
            if DEBUG:
                print('STOPPING at pair', k, ':', g1, '+', g2, 'metric =', pair_metric)
            break

        if DEBUG:
            print('Testing pair', k, ':', g1, '+', g2, 'metric =', pair_metric)

        num_swap_temp = 0

        central_pair, swap_pair = CA_SP_det(g1, g2, pref_pair)

        indx_first=-1
        indx_second=-1
        indx_first_swap=-1
        indx_second_swap=-1

        for num,item in enumerate(elements):
            if(item == central_pair[0]):
                indx_first=num
            if(item == central_pair[1]):
                indx_second=num
        for num,item in enumerate(elements):
            if(item == swap_pair[0]):
                indx_first_swap=num
            if(item == swap_pair[1]):
                indx_second_swap=num

        if(indx_first == 0):
            low_range=0
            high_range=num_each[indx_first]
        else:
            low_range=sum(num_each[a] for a in range(indx_first))
            high_range=sum(num_each[a] for a in range(indx_first+1))

        if DEBUG:
            print('low_range,high_range={},{}'.format(low_range,high_range))

        swap_found_this_pair = False

        for a in range(low_range,high_range):
            exist=True
            if (len(forbidden_vals[indx_first]) == 0):
                exist=False
            else:
                for item in forbidden_vals[indx_first]:
                    if(item == a):
                        exist=True
                        break
                    else:
                        exist=False
            if(exist == True):
                continue
            elif(exist==False):
                avl_space=num_each[indx_first]-counts[indx_first]
                total_CA1_selected += 1
                if (avl_space >= 1):
                    nnlist1=nn_dict[a]
                    forbidden_list.append(a)
                    forbidden_vals[indx_first][int(counts[indx_first])]=a
                    counts[indx_first]+=1
                else:
                    terminated_due_to_full += 1
                    if DEBUG:
                        print("FULL CONDITION", 'full', len(forbidden_vals[indx_first]))
                    break
            if DEBUG and total_CA1_selected % 100 == 0:
                print(total_CA1_selected, successful_swaps, failed_swaps)
            b=True
            e=True
            while b == True or e == True:
                c=random.randint(num_each[indx_second-1],(num_each[indx_second-1]+num_each[indx_second]-1))
                for d in nnlist1:
                    if(d==c):
                        b=True
                        break
                    else:
                        b=False
                if(len(forbidden_vals[indx_second]) == 0):
                    e=False
                else:
                    for item in forbidden_vals[indx_second]:
                        if(item == c):
                            e=True
                            break
                        else:
                            e=False
            avl_space=num_each[indx_second]-counts[indx_second]
            if(avl_space >= 1):
                forbidden_list.append(c)
                forbidden_vals[indx_second][int (counts[indx_second])]=c
                counts[indx_second]+=1
                nnlist2=nn_dict[c]
            else:
                terminated_due_to_full += 1
                break
            num_each_max=np.full(len(elements),0,dtype=int)
            num_each_counter=0
            for x1 in range(len(elements)):
                num_each_max[x1]=num_each[x1]+num_each_counter
                num_each_counter=num_each_max[x1]
            present1=False
            present2=True
            swapable_first=[]
            for f in nnlist1:
                for g in forbidden_list:
                    if(f==g):
                        present1=True
                        break
                    else:
                        present1=False
                if(indx_first_swap == 0):
                    if(f < num_each[indx_first_swap]):
                        present2=False
                    else:
                        present2=True
                else:
                    if(f >= num_each_max[indx_first_swap-1] and f < num_each_max[indx_first_swap]):
                        present2=False
                    else:
                        present2=True
                if present1 == False and present2 == False:
                    swapable_first.append(f)

            present3=False
            present4=True
            swapable_second=[]
            for m in nnlist2:
                for n in forbidden_list:
                    if(m==n):
                        present3=True
                        break
                    else:
                        present3=False
                if(indx_second_swap == 0):
                    if(m < num_each[indx_second_swap]):
                        present4=False
                    else:
                        present4=True
                else:
                    if(m >= num_each_max[indx_second_swap-1] and  m < num_each_max[indx_second_swap]):
                        present4=False
                    else:
                        present4=True
                if present3 == False and present4 == False:
                    swapable_second.append(m)

            id_central_pref=-1
            for p in pref_pair:
                for q,r in enumerate(central_pair):
                    if(p==r):
                        id_central_pref=q
                        break

            if(id_central_pref == 0):
                for item in nnlist1:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_second_swap):
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break
                for num,item1 in enumerate(swapable_second):
                    nnlist=nn_dict[item1]
                    for item2 in nnlist:
                        elem_indx=func_all.indx_find(num_each,item2,elements)
                        if(elem_indx == indx_first):
                            avl_space=num_each[elem_indx]-counts[elem_indx]
                            if(avl_space >= 1):
                                forbidden_vals[elem_indx][int (counts[elem_indx])]=item1
                                counts[elem_indx]+=1
                                swapable_second.pop(num)
                                forbidden_list.append(item2)
                                break
                            else:
                                break
            elif(id_central_pref == 1):
                for item in nnlist2:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_first_swap):
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break
                for num,item1 in enumerate(swapable_first):
                    nnlist=nn_dict[item1]
                    for item2 in nnlist:
                        elem_indx=func_all.indx_find(num_each,item2,elements)
                        if(elem_indx == indx_second):
                            avl_space=num_each[elem_indx]-counts[elem_indx]
                            if(avl_space >= 1):
                                forbidden_vals[elem_indx][int (counts[elem_indx])]=item1
                                counts[elem_indx]+=1
                                swapable_first.pop(num)
                                forbidden_list.append(item2)
                                break
                            else:
                                break
            else:
                raise Exception('id_central_pair could not be determined')

            if DEBUG:
                print('a,swapable_first,swapable_second,counts[indx_first_swap],counts[indx_second_swap]={},{},{},{},{}'.format(
                    a, swapable_first, swapable_second, counts[indx_first_swap], counts[indx_second_swap]))

            if(len(swapable_first) >= len(swapable_second) and len(swapable_first) != 0 and len(swapable_second) != 0):
                successful_swaps += 1
                max_swap_pos=-1
                temp_arr=np.zeros(3)
                temp_arr[0]=num_each[indx_first_swap]-counts[indx_first_swap]
                temp_arr[1]=num_each[indx_second_swap]-counts[indx_second_swap]
                temp_arr[2]=len(swapable_second)
                max_swap_pos=int (np.min(temp_arr))
                for num in range(max_swap_pos):
                    forbidden_vals[indx_first_swap][int (counts[indx_first_swap])]=swapable_first[num]
                    counts[indx_first_swap]+=1
                    forbidden_list.append(swapable_first[num])
                    forbidden_vals[indx_second_swap][int (counts[indx_second_swap])]=swapable_second[num]
                    counts[indx_second_swap]+=1
                    forbidden_list.append(swapable_second[num])
                    swap_dict[swapable_first[num]]=swapable_second[num]
                    indx_red1=func_all.pair_index(pairs,central_pair[0],swap_pair[0])
                    indx_red2=func_all.pair_index(pairs,central_pair[1],swap_pair[1])
                    indx_inc1=func_all.pair_index(pairs,central_pair[0],swap_pair[1])
                    indx_inc2=func_all.pair_index(pairs,central_pair[1],swap_pair[0])
                    num_swap+=1
                    num_swap_temp+=1
                swap_found_this_pair = True
                # break

            elif(len(swapable_first) < len(swapable_second) and len(swapable_first) != 0 and len(swapable_second) != 0):
                successful_swaps += 1
                max_swap_pos=-1
                temp_arr=np.zeros(3)
                temp_arr[0]=num_each[indx_first_swap]-counts[indx_first_swap]
                temp_arr[1]=num_each[indx_second_swap]-counts[indx_second_swap]
                temp_arr[2]=len(swapable_first)
                max_swap_pos=int (np.min(temp_arr))
                for num in range(max_swap_pos):
                    forbidden_vals[indx_first_swap][int (counts[indx_first_swap])]=swapable_first[num]
                    counts[indx_first_swap]+=1
                    forbidden_list.append(swapable_first[num])
                    forbidden_vals[indx_second_swap][int (counts[indx_second_swap])]=swapable_second[num]
                    counts[indx_second_swap]+=1
                    forbidden_list.append(swapable_second[num])
                    swap_dict[swapable_first[num]]=swapable_second[num]
                    indx_red1=func_all.pair_index(pairs,central_pair[0],swap_pair[0])
                    indx_red2=func_all.pair_index(pairs,central_pair[1],swap_pair[1])
                    indx_inc1=func_all.pair_index(pairs,central_pair[0],swap_pair[1])
                    indx_inc2=func_all.pair_index(pairs,central_pair[1],swap_pair[0])
                    num_swap+=1
                    num_swap_temp+=1
                swap_found_this_pair = True
                # break

            elif (len(swapable_first) == 0 or len(swapable_second) == 0):
                continue

        if swap_found_this_pair:
            pairs_with_swap += 1
            num_bond_mod[indx_red1]=num_bond_mod[indx_red1]-num_swap_temp
            num_bond_mod[indx_red2]=num_bond_mod[indx_red2]-num_swap_temp
            num_bond_mod[indx_inc1]=num_bond_mod[indx_inc1]+num_swap_temp
            num_bond_mod[indx_inc2]=num_bond_mod[indx_inc2]+num_swap_temp
        else:
            pairs_without_swap += 1
            failed_swaps += 1
            if DEBUG:
                print("NO SWAP FOUND FOR PAIR:", g1, g2)

    count_bonds_mod=[]
    for num,item in enumerate(count_bonds):
        count_bonds_mod.append(item+num_bond_mod[num])

    print("\n" + "="*60)
    print("           swap_pos() Performance Summary")
    print("="*60)
    print(f"Total CA1 selections        : {total_CA1_selected}")
    print(f"Successful swaps            : {successful_swaps}")
    print(f"Pairs tested (matched, sorted) : {total_pairs_considered} / {len(group1)}")
    print(f"Pairs that produced a swap   : {pairs_with_swap}")
    print(f"Pairs with no swap           : {pairs_without_swap}")
    print(f"Terminated because full      : {terminated_due_to_full}")
    print("="*60)

    return num_swap,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2