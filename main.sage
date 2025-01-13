print("using RSK correspondence to calculate Lusztig's a function")


#######################################################################################

## --- python packages for plotting the results

#######################################################################################


from collections import Counter
from collections import defaultdict
from tabulate import tabulate
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt

import itertools



#######################################################################################

## --- sagemath uses another set of generator T_i = qH_i for the Iwahori-Hecke algebra if we set the scalars as follows

#######################################################################################

n = 4

R.<q> = LaurentPolynomialRing(QQ)
scalar = q

H = IwahoriHeckeAlgebra(f'A{n-1}', -scalar^2,1)
T = H.T()


tgen = T.algebra_generators()
H_simple = [t/scalar for t in tgen]


def standard_basis_element_from_word(list_of_indices):
    out = 1
    for i in list_of_indices:
        out *= H_simple[i-1]  
    return out


#######################################################################################

## --- Writing the Hecke Welements in my normalization

#######################################################################################

def ordering_Sn_elements(n):
    ##### reduced words in the simple reflections, in the order of the T basis Welements in sage
    perms = list(WeylGroup(['A',n-1],prefix='s'))
    perms = [perm.reduced_word() for perm in perms]
    perms = sorted(perms,key=lambda x: (len(x),x))
    return perms



## just a formatting function
def print_in_basis(vec, basis_words, basis_name: str):   
    display = ''
    for coeff, word in zip(vec,basis_words):
        if coeff!=0: 
            display += f"({coeff})*{basis_name}{word} + " if word else f"{coeff} + "
    print(display[:-2])




def hecke_in_my_normalization(t,detailed=False):###### converts the T basis to the H basis
    #basis = t.parent().gens()
    coeffs = t.to_vector()

    ##### now need to figure out which q power to multiply by, namely the length of the corresponding permutation

    gens = list(t.parent().algebra_generators())
    n = len(gens) +1 ##### size of Sn

    perms = ordering_Sn_elements(n)
    out = [0]*len(perms)
    for i,(coeff,word) in enumerate(zip(coeffs,perms)):
        if coeff!=0: out[i] = coeff*(scalar**len(word))

    if detailed: 
        #### now print the result in terms of the simple reflections
        print_in_basis(out, perms, 'H')

    return out


print('\n\n\nthe quadratic relation in my normalization: H[1]^2 == ')
t = hecke_in_my_normalization(H_simple[0]^2,True) ### normalization check


#######################################################################################

## --- now let us check the normalization of the KL basis Welements

#######################################################################################


KL_basis = H.C()
KL_basis_list = list(KL_basis.basis())



print('\n\n\nthe KL basis element of H[1] in my normalization: ')
t = hecke_in_my_normalization(T(KL_basis[1]),True) 

print('\n\n\nthe KL basis element of H[1,2,1] in my normalization: ')
t = hecke_in_my_normalization( T(KL_basis_list[-1]),True)




#######################################################################################

## --- compute the a function for the symmetric group S_n a la lusztig's definition

#######################################################################################



#kl_polynom_matrix = matrix([hecke_in_my_normalization(T(x)) for x in KL_basis_list])
#### KL basis in terms of the standard H basis

W = WeylGroup(f"A{n-1}")
Welements= W.list()
Welements = list(sorted(Welements,key=lambda x: (x.length(),x.reduced_word())))

print(f'\n\n\nelements of S_{n} as permutations:\n', [x.to_permutation() for x in Welements])

e = Welements[0]
s1=Welements[1]
s2=Welements[2]

def h(x,y,z,detailed=False):

    Cx = KL_basis_list[Welements.index(x)]
    Cy = KL_basis_list[Welements.index(y)]
    product = KL_basis( Cx * Cy ).to_vector()
    if detailed: 
        print_in_basis(product, W, 'C')
    return product[Welements.index(z)]



def a_function_lusztig(z):
    W = z.parent()
    aux = [ h(x,y,z) for x in W for y in W]
    #print('z,aux:',z.to_permutation(),aux)
    aux = [R(f).dict() for f in aux]
    aux = [power for dic in aux for power in dic.keys()]
    #print('z,aux:',z.to_permutation(),aux)

    min_power = min(aux) if aux else 0  # minimum power of q in the coefficients
    return -min_power if min_power<0 else 0




#######################################################################################

## --- now let us compute the a function for the symmetric group S_n via RSK

#######################################################################################

def n_choose2(n):
    return n*(n-1)//2

def a_function(w:list):
    rsk=RSK(w)
    Q=rsk[1].conjugate()
    return sum([n_choose2(len(z))for z in Q])


def get_multiplicities(lst):
    return Counter(lst)


#### check that this agrees with the previous a function via Lusztigs definition, but h is slow for n>=4

def check_a_func_agree():
    print(f'\n\n\na function values via the two definitions for elements of S_{n}:\n')
    data = []
    for z in Welements: 
        zz = z.to_permutation()
        a1 = a_function_lusztig(z)
        a2 = a_function(zz)
        print('permutation, a_function_lusztig, a_function_RSK:',zz,a1,a2)
        data.append((zz,a1,a2))

        if a1 != a2: raise ValueError('the two definitions of the a function do not agree')
    #headers = ['permutation', 'a_function_lusztig', 'a_function_RSK']
    #print(tabulate(data, headers=headers, tablefmt='grid'))

#check_a_func_agree()

#######################################################################################

## --- a function values for all Welements in S_n, via RSK

#######################################################################################

def a_function_table(n):
        
    W = WeylGroup(f"A{n-1}")
    # [s1,s2,s3] = W.simple_reflections()
    KL = KazhdanLusztigPolynomial(W,q)
    Welements= W.list()
    KLpolys = [KL.P(Welements[0],s) for s in Welements]
    sn_elements = [w.to_permutation() for w in Welements]


    out= [a_function(w) for w in sn_elements]
    lengths = [Permutation(list(w)).number_of_inversions() for w in sn_elements]
    count = get_multiplicities(out)
    #print('values and their multiplicities:')
    #print(count)
    #print('total number of values:', len(count))
    #print(f'number of partitions of {n}:',number_of_partitions(n) )
    zipped_list = list(zip(sn_elements, KLpolys, out, lengths))
    sorted_zipped_list = sorted(zipped_list, key=lambda x: (x[-2], x[-1]))
    
    D = []
    # Prepare data for tabulate
    table_data = []
    for w,P, a, l in sorted_zipped_list:
        aux = l - 2 * P.degree()
        row = [w, P, a,  aux,l]
        if aux==a: 
            D.append((w,P,a,l)) 
            row = [f"\033[1;31m{cell}\033[0m" for cell in row]  # Red color
        table_data.append(row)

    headers = ['w', 'P(e,w)', 'a(w)', 'l(w) - 2 deg(P)', 'l(w)']

    print(f'\n\n\ntable of a-function values for S_{n}')
    print(tabulate(table_data, headers=headers, tablefmt='grid'))

    return sorted_zipped_list,D

data, D = sorted_result = a_function_table(n)



#######################################################################################

## --- distinguished involution as max min coset reps?

#######################################################################################




def parabolics(n):
    original_set = [i for i in range(1,n)]

    # Generate all subsets
    all_subsets = []
    for r in range(len(original_set) + 1):
        subsets = list(itertools.combinations(original_set, r))
        all_subsets.extend(subsets)
    return all_subsets



def max_min_coset_reps(n):
    W = WeylGroup(f"A{n-1}")
    Welements= W.list()
    parabolics_list = parabolics(n)
    data = [[w.coset_representative(p) for w in Welements] for p in parabolics_list]
    data = [sorted(x,key=lambda x: (x.length(),x.reduced_word())) for x in data]
    data = [x[-1].to_permutation() for x in data]
    return data

print(max_min_coset_reps(3))#### this is not right