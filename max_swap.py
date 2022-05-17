import numpy as np
import random
import func_all

def swap_pos(pref_pair_req,pairs,elements,num_each,nn_dict,count_bonds):
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
        print('Issue with pref_pair list')
        exit(1)
    #intialising a 2D array, which would store the forbidden indices
    rows=len(elements)
    #keeping a chance for the cases, where number of elements might be different!
    cols=[]
    for item in num_each:
        cols.append(item)
    forbidden_vals=[[-1 for a in range(cols[b])] for b in range(rows)]
    counts=np.zeros(len(elements))
    forbidden_list=[]
    num_swap=0
    num_bond_mod=np.zeros(len(pairs)) #array to store the number of bons increase/decrease
    group1=[]
    group2=[]
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
                    group1.append(item)
            if(temp1 == pref_pair[1] or temp2 == pref_pair[1]): #pref_pair[1] --> group2
                if((temp1 == pref_pair[0] and temp2 == pref_pair[1]) or (temp1==pref_pair[1] and temp2==pref_pair[0])):
                    pass #pref_pair is not added in the group1
                else:
                    group2.append(item)
        if(pref_pair[0] == pref_pair[1]):
            temp1=item.split('-')[0]
            temp2=item.split('-')[1]
            if(temp1==pref_pair[0] or temp2==pref_pair[1]):
                #print(temp1,temp2)
                if(temp1==pref_pair[0] and temp2==pref_pair[1]):
                    pass
                elif(count_like%2 != 0):
                    group1.append(item)
                    count_like+=1
                elif(count_like%2 == 0):
                    group2.append(item)
                    count_like+=1
    
    print (group1,group2)
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
            #for num in range(temp):
            #    group1.pop(num)
        elif(len(group1) < len(group2)):
            temp=len(group2)-len(group1)
            #temp should always be equal to 1
            #as to make any number even from odd only one should be removed
            if (temp != 1):
                print ('Issue with group1 and group2 generation-2')
                exit(1)
            else:
                pass
            num=random.randint(0,(len(group2)-1))
            group2.pop(num)

            #for num in range(temp):
            #    group2.pop(num)

    #Initialisation of swap_dict: key is the first atom and val is the second atom, involved in the swap.
    swap_dict={}
    #TODO at this moment, the pref_pair and associated change are determined
    #on the basis of sequential access to group-1. We need to take bond from
    #group1, such that prespecified pref_bond instances would increase, and 
    #and unpref_bond (which needs to be specied) instances would decrease.
    for x in group1:
        num_swap_temp=0
        swap_pair=[]
        central_pair=[]
        first_x=x.split('-')[0]
        second_x=x.split('-')[1]
        if(first_x == pref_pair[0] or first_x == pref_pair[1]):
            unpref_at1=second_x
            pref_at1=first_x
        elif(second_x == pref_pair[0] or second_x == pref_pair[1]):
            pref_at1=second_x
            unpref_at1=first_x
        else:
            print('Issue with the unpref_at1 determination')
        #print(first_x,second_x,pref_pair[0],pref_pair[1],unpref_at1,pref_at1)
        for num_y,y in enumerate(group2):
            first_y=y.split('-')[0]
            second_y=y.split('-')[1]
            if(first_y == pref_pair[0] or first_y == pref_pair[1]):
                pref_at2=first_y
                unpref_at2=second_y
            elif(second_y == pref_pair[0] or second_y == pref_pair[1]):
                pref_at2=second_y
                unpref_at2=first_y
            else:
                print('Issue with pref_at2 determination')
            
            if(pref_pair[0] != pref_pair[1]): #and unpref_at1 != unpref_at2): #unlike atoms pref
                swap_pair.append(unpref_at1)
                swap_pair.append(pref_at2)
                central_pair.append(pref_at1)
                central_pair.append(unpref_at2)
                group2.pop(num_y)
                break
            elif(pref_pair[0] == pref_pair[1]): #and unpref_at1 != unpref_at2):
                swap_pair.append(unpref_at1)
                swap_pair.append(pref_at1)
                central_pair.append(pref_at1)
                central_pair.append(unpref_at2)
                group2.pop(num_y)
                break
        print('swap_pair,central_pair={}{}'.format(swap_pair,central_pair))

        #finding the maximum like atom swap possible
        #determination of indices of pref_pair, as defined in the element list
        #swap_pair=['Mn','Mn']
        #central_pair=['Ni','Ni']
        #pref_pair=['Mn','Ni'] 
        indx_first=-1 #index of first element of the central_pair
        indx_second=-1 #index of second element of the central_pair
        indx_first_swap=-1 #index of first element of the swap_pair
        indx_second_swap=-1 #index of second element of the swap_pair

        for num,item in enumerate(elements):
            if(item == central_pair[0]):
                indx_first=num
            if(item == central_pair[1]):
                indx_second=num
        print('indx_first,indx_second={},{}'.format(indx_first,indx_second))
        for num,item in enumerate(elements):
            if(item == swap_pair[0]):
                indx_first_swap=num
            if(item == swap_pair[1]):
                indx_second_swap=num
        print('indx_first_swap,indx_second_swap={},{}'.format(indx_first_swap,indx_second_swap))

        #intialising a 2D array, which would store the forbidden indices
        #rows=len(elements)
        #keeping a chance for the cases, where number of elements might be different!
        #cols=[]
        #for item in num_each:
        #    cols.append(item)
        #forbidden_vals=[[-1 for a in range(cols[b])] for b in range(rows)]
        #counts=np.zeros(len(elements))
        #forbidden_list=[]
        #num_swap=0
        if(indx_first == 0):
            low_range=0
            high_range=num_each[indx_first]
        else:
            low_range=sum(num_each[a] for a in range(indx_first))
            high_range=sum(num_each[a] for a in range(indx_first+1))

        #print(low_range,high_range)
        for a in range(low_range,high_range): #num_each[0]
            exist=True  #binary variable depicting whether the central_pair[0] is available.
            #exist=True; If a is in the forbidden_list, False; otherwise.
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
            elif(exist==False): #if a is not in the forbidden_list already
                #ensuring that addition is not allowed beyond max limit
                avl_space=num_each[indx_first]-counts[indx_first]
                if (avl_space >= 1):
                    nnlist1=nn_dict[a]
                    forbidden_list.append(a)
                    forbidden_vals[indx_first][int(counts[indx_first])]=a
                    counts[indx_first]+=1
                else:
                    #print('full',len(forbidden_vals[indx_first]))
                    break
            b=True #boolean parameter depicting whether randomly chosen central_pair[1] ele is already in nnlist1 
            e=True #boolean parameter depicting whether randomly chosen cental_pair[1] ele is in forbidden_list 
            while b == True or e == True:
                #central_pair[1] is being randomly chosen.
                c=random.randint(num_each[indx_second-1],(num_each[indx_second-1]+num_each[indx_second]-1))
                for d in nnlist1:
                    if(d==c):
                        b=True #c is in nnlist1, so other c needs to be chosen.
                        break
                    else:
                        b=False
                #checking in the forbidden list
                if(len(forbidden_vals[indx_second]) == 0):
                    e=False
                else:
                    for item in forbidden_vals[indx_second]:
                        if(item == c):
                            e=True
                            break
                        else:
                            e=False
            #ensuring the addition is not allowed beyond max limit
            avl_space=num_each[indx_second]-counts[indx_second]
            if(avl_space >= 1):
                forbidden_list.append(c)
                forbidden_vals[indx_second][int (counts[indx_second])]=c
                counts[indx_second]+=1
                nnlist2=nn_dict[c]
            else:
                break
            num_each_max=np.full(len(elements),0,dtype=int) #array to store the max index for each ele
            num_each_counter=0
            for x in range(len(elements)):
                    num_each_max[x]=num_each[x]+num_each_counter
                    num_each_counter=num_each_max[x]
            present1=False #boolean para depicting if the nnlist1 element is in forbidden_list
            present2=True  #boolean parameter depicting if ele in nnlist is swap_pair[0]
            swapable_first=[] #list to store the ele in nnlist1, which may be swapped.
            #In this part, elements in nnlist1 are determined, which may be swapped with ele
            #of nnlist2.
            for f in nnlist1:
                for g in forbidden_list:
                    if(f==g):
                        present1=True
                        break
                    else:
                        present1=False #element is not in FL
                if(indx_first_swap == 0):
                    if(f < num_each[indx_first_swap]):
                        present2=False
                    else:
                        present2=True
                else:
                    if(f >= num_each_max[indx_first_swap-1] and f < num_each_max[indx_first_swap]):
                        #print ('checking swapable_first={},{},{}'.format(f,num_each_max[indx_first_swap-1],\
                        #        num_each_max[indx_first_swap]))
                        present2=False #ele in nnlist1 is swap_pair[0]
                    else:
                        present2=True

                if present1 == False and present2 == False:
                    swapable_first.append(f)
                    #counts[indx_first_swap]+=1
                
            present3=False #boolean para depicting if nnlist2 element is in the FL
            present4=True #boolean para depicting, if ele in nnlist2 is swap_pair[1] 
            swapable_second=[] #list to store the ele in nnlist2, which may be swapped
            #with the ele in nnlist1
            for m in nnlist2:
                for n in forbidden_list:
                    if(m==n):
                        present3=True
                        break
                    else:
                        present3=False #ele is not in FL
                if(indx_second_swap == 0):
                    if(m < num_each[indx_second_swap]):
                        present4=False
                    else:
                        present4=True
                else:
                    if(m >= num_each_max[indx_second_swap-1] and  m < num_each_max[indx_second_swap]):
                        #print('m value for swapable_second={}'.format(m))
                        present4=False
                    else:
                        present4=True
                if present3 == False and present4 == False:
                    swapable_second.append(m)
            #Depending upon the existance of the central atom in the pref_pair, algorithm
            #for the swap need to be devised. 
            id_central_pref=-1
            for p in pref_pair:
                for q,r in enumerate(central_pair):
                    if(p==r):
                        id_central_pref=q
                        break
            #print ('ide_central_pref={}'.format(id_central_pref)) 
            if(id_central_pref == 0):
                #to ensure that pref_bond is not removed, nnlist1 is checked for 
                #swap_pair[1] atom and if found, it is frozen for the swap and 
                #interred into the FL.
                for item in nnlist1:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_second_swap): #if,desired atom is already in the nnlist1, freeze it!
                        #ensuring that addtion beyond max limit is not allowed.
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break
                #Similarly, for nnlist2, if swap_pair[0] is found, it
                #needs to be forzen and entered into FL.
                for item in nnlist2:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_first_swap): #if,desired atom is already in the nnlist1, freeze it!
                        #ensuring that addtion beyond max limit is not allowed.
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break

                #print(a,len(swapable_second))
                for num,item1 in enumerate(swapable_second):
                    nnlist=nn_dict[item1]
                    for item2 in nnlist:
                        elem_indx=func_all.indx_find(num_each,item2,elements)
                        if(elem_indx == indx_first): #if the NN of the element in swapable_second is 
                        #already central_pair[0], then such element is already forming pref_bond, so it 
                        #should be frozen.
                            #ensuring that addtion beyond max limit is not allowed
                            avl_space=num_each[elem_indx]-counts[elem_indx]
                            if(avl_space >= 1):
                                forbidden_vals[elem_indx][int (counts[elem_indx])]=item1
                                counts[elem_indx]+=1
                                swapable_second.pop(num)
                                forbidden_list.append(item2)
                                break
                            else:
                                break

            if(id_central_pref == 1): #if central_pair[1] is in the pref_pair
                #freezing swap_pair[0] atoms in nnlist2 
                for item in nnlist2:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_first_swap):
                        #ensuring that addtion beyond max limit is not allowed.
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break
                #freezing swap_pair[1] atoms in nnlist1 
                for item in nnlist1:
                    indx=func_all.indx_find(num_each,item,elements)
                    if(indx == indx_second_swap):
                        #ensuring that addtion beyond max limit is not allowed.
                        avl_space=num_each[indx]-counts[indx]
                        if(avl_space >= 1):
                            forbidden_vals[indx][int (counts[indx])]=item
                            counts[indx]+=1
                            forbidden_list.append(item)
                        else:
                            break

                #print(a,len(swapable_first))
                for num,item1 in enumerate(swapable_first):
                    nnlist=nn_dict[item1]
                    for item2 in nnlist:
                        elem_indx=func_all.indx_find(num_each,item2,elements)
                        if(elem_indx == indx_second):
                            #ensuring that addtion beyond max limit is not allowed
                            avl_space=num_each[elem_indx]-counts[elem_indx]
                            if(avl_space >= 1):
                                forbidden_vals[elem_indx][int (counts[elem_indx])]=item1
                                counts[elem_indx]+=1
                                swapable_first.pop(num)
                                forbidden_list.append(item2)
                                break
                            else:
                                break

            print(a,len(swapable_first),len(swapable_second),counts[indx_first_swap],counts[indx_second_swap])
            if(len(swapable_first) >= len(swapable_second) and len(swapable_first) != 0 and len(swapable_second) != 0):
                #print(counts[indx_first_swap],counts[indx_second_swap])
                #ensuring that addition beyond max limit is not allowed
                max_swap_pos=-1
                temp_arr=np.zeros(3) #this array will store the 3 values, of ehich min would be taken
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
                    #SWAP: Here the identities of the atoms, which may be swapped is being stored
                    #we are not carrying out swap here, as we would like to control the number of
                    #swaps and hence, the value of the order parameter. we would generate a dict,
                    #whose, key would be the ID of the first atom and value would be second atom
                    #involved in the swap. swap_dict{} would be initialised before the outermost 
                    #loop.
                    swap_dict[swapable_first[num]]=swapable_second[num]
                    #In this part, the pair indices, which are being reduced and increased are 
                    #being determined, with an ASSUMPTION, that number of swaps with composition 
                    #constraint leads to the situation, in which as we are swapping the atoms
                    #and updating the forbidden_list, it get ppopulated fully with only one pair
                    #of bond swap, we may not need any other bond pairs from group1 and group2.
                    indx_red1=func_all.pair_index(pairs,central_pair[0],swap_pair[0])
                    indx_red2=func_all.pair_index(pairs,central_pair[1],swap_pair[1])
                    indx_inc1=func_all.pair_index(pairs,central_pair[0],swap_pair[1])
                    indx_inc2=func_all.pair_index(pairs,central_pair[1],swap_pair[0])

                    #print('one',rand_pos[swapable_first[num]],rand_pos[swapable_second[num]])
                    #temp1=rand_pos[swapable_first[num]]
                    #temp2=rand_pos[swapable_second[num]]
                    #rand_pos[swapable_first[num]]=temp2
                    #rand_pos[swapable_second[num]]=temp1
                    #print('second',rand_pos[swapable_first[num]],rand_pos[swapable_second[num]])
                    num_swap+=1
                    num_swap_temp+=1
            elif(len(swapable_first) < len(swapable_second) and len(swapable_first) != 0 and len(swapable_second) != 0):
                #ensuring that addition beyond max limit is not allowed
                max_swap_pos=-1
                temp_arr=np.zeros(3) #this array will store the 3 values, of ehich min would be taken
                temp_arr[0]=num_each[indx_first_swap]-counts[indx_first_swap]
                temp_arr[1]=num_each[indx_second_swap]-counts[indx_second_swap]
                temp_arr[2]=len(swapable_first)
                max_swap_pos=int (np.min(temp_arr))
            
                for num in range(max_swap_pos):
                    #print(indx_first_swap,counts[indx_first_swap],num,len(swapable_first))
                    #id of swap, swapable_first[num] and swapable_second[num 
                    #print(swapable_first[num],swapable_second[num])
                    forbidden_vals[indx_first_swap][int (counts[indx_first_swap])]=swapable_first[num]
                    counts[indx_first_swap]+=1
                    forbidden_list.append(swapable_first[num])
                    forbidden_vals[indx_second_swap][int (counts[indx_second_swap])]=swapable_second[num]
                    counts[indx_second_swap]+=1
                    forbidden_list.append(swapable_second[num])
                    #SWAP: Here the identities of the atoms, which may be swapped is being stored
                    #we are not carrying out swap here, as we would like to control the number of
                    #swaps and hence, the value of the order parameter. we would generate a dict,
                    #whose, key would be the ID of the first atom and value would be second atom
                    #involved in the swap. swap_dict{} would be initialised before the outermost 
                    #loop.
                    swap_dict[swapable_first[num]]=swapable_second[num]
                    #In this part, the pair indices, which are being reduced and increased are 
                    #being determined, with an ASSUMPTION, that number of swaps with composition 
                    #constraint leads to the situation, in which as we are swapping the atoms
                    #and updating the forbidden_list, it get ppopulated fully with only one pair
                    #of bond swap, we may not need any other bond pairs from group1 and group2.
                    indx_red1=func_all.pair_index(pairs,central_pair[0],swap_pair[0])
                    indx_red2=func_all.pair_index(pairs,central_pair[1],swap_pair[1])
                    indx_inc1=func_all.pair_index(pairs,central_pair[0],swap_pair[1])
                    indx_inc2=func_all.pair_index(pairs,central_pair[1],swap_pair[0])

                    #print(rand_pos[swapable_first[num]],rand_pos[swapable_second[num]])
                    #temp1=rand_pos[swapable_first[num]]
                    #temp2=rand_pos[swapable_second[num]]
                    #rand_pos[swapable_first[num]]=temp2
                    #rand_pos[swapable_second[num]]=temp1
                    #print(rand_pos[swapable_first[num]],rand_pos[swapable_second[num]])
                    num_swap+=1
                    num_swap_temp+=1

            #print(a,counts[indx_first_swap],counts[indx_second_swap],num_swap)
        #print(central_pair,swap_pair,num_swap_temp)
        #indx_red1=func_all.pair_index(pairs,central_pair[0],swap_pair[0])
        #indx_red2=func_all.pair_index(pairs,central_pair[1],swap_pair[1])
        #indx_inc1=func_all.pair_index(pairs,central_pair[0],swap_pair[1])
        #indx_inc2=func_all.pair_index(pairs,central_pair[1],swap_pair[0])
        print('pairs={}'.format(pairs))
        print('indx_red1,indx_red2,indx_inc1,indx_inc2={},{},{},{}'.format(indx_red1,indx_red2,indx_inc1,indx_inc2))
        #print(pairs)
        #updating the reduced and increased bonds
        num_bond_mod[indx_red1]=num_bond_mod[indx_red1]-num_swap_temp
        num_bond_mod[indx_red2]=num_bond_mod[indx_red2]-num_swap_temp
        num_bond_mod[indx_inc1]=num_bond_mod[indx_inc1]+num_swap_temp
        num_bond_mod[indx_inc2]=num_bond_mod[indx_inc2]+num_swap_temp
    #print(num_swap)
    #print(num_bond_mod)
    count_bonds_mod=[]
    for num,item in enumerate(count_bonds):
        count_bonds_mod.append(item+num_bond_mod[num])
    return num_swap,swap_dict,indx_red1,indx_red2,indx_inc1,indx_inc2