print("using RSK correspondence to calculate Lusztig's a function")

reset()
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

## --- Global variables: weyl group, hecke algebra, KL basis elements

#######################################################################################

n = 4

R.<v> = LaurentPolynomialRing(QQ)


###### Weyl group
W = WeylGroup(f"A{n-1}")
Welements= W.list()
#Welements = list(sorted(Welements,key=lambda x: (x.length(),x.reduced_word())))
#this is the order of the basis elements if below we use IwahoriHeckeAlgebra('A3') instead of W, I hate this

print(f'\n\n\nelements of S_{n} as permutations:\n', [x.to_permutation() for x in Welements])

e = Welements[0]
s1=Welements[1]
s2=Welements[2]

##### hecke algebra, sagemath uses another set of generator T_i = qH_i for the Iwahori-Hecke algebra if we set the scalars as follows, where H_i is our usual generator
scalar = v

H = IwahoriHeckeAlgebra(W, -scalar^2,1)
T = H.T()


tgen = T.algebra_generators()
H_simple = [t/scalar for t in tgen]


#######################################################################################

## --- Writing the Hecke elements in my normalization

#######################################################################################


def standard_basis_element_from_word(list_of_indices):
    out = 1
    for i in list_of_indices:
        out *= H_simple[i-1]  
    return out


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

    ##### now need to figure out which v power to multiply by, namely the length of the corresponding permutation

    gens = list(t.parent().algebra_generators())
    #n = len(gens) +1 ##### size of Sn
    #perms = ordering_Sn_elements(n)
    perms = [x.reduced_word() for x in Welements]
    out = [0]*len(perms)
    for i,(coeff,word) in enumerate(zip(coeffs,perms)):
        if coeff!=0: out[i] = coeff*(scalar**len(word))

    if detailed: 
        #### now print the result in terms of the simple reflections
        print_in_basis(out, perms, 'H')

    return out


print('\n\n\nthe quadratic relation in my normalization: H[1]^2 == ')
t = hecke_in_my_normalization(H_simple[0]^2,True) ### normalization check

#### stanard basis elements in H(S3)
h1 = standard_basis_element_from_word([1])
h2 = standard_basis_element_from_word([2])
h12 = standard_basis_element_from_word([1,2])
h21 = standard_basis_element_from_word([2,1])
h121 = standard_basis_element_from_word([1,2,1])

##### KL basis elements in H(S3)
c1 = standard_basis_element_from_word([1])+v
c2 = standard_basis_element_from_word([2])+v
c12 = c1*c2
c21 = c2*c1
c121 = c12*c1 - c1


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

    min_power = min(aux) if aux else 0  # minimum power of v in the coefficients
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
    print(f'\n\n\na function values via the two definitions for W_elements of S_{n}:\n')
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
    KL = KazhdanLusztigPolynomial(W,v)

    sn_elements= W.list()
    KLpolys = [KL.P(Welements[0],s) for s in Welements]

    out= [a_function(w.to_permutation()) for w in sn_elements]
    
    count = get_multiplicities(out)
    #print('values and their multiplicities:')
    #print(count)
    #print('total number of values:', len(count))
    #print(f'number of partitions of {n}:',number_of_partitions(n) )
    zipped_list = list(zip(sn_elements, KLpolys, out))
    sorted_zipped_list = sorted(zipped_list, key=lambda x: (x[-1], x[0].length()))
    
    out = []
    D = []
    # Prepare data for tabulate
    table_data = []
    for w,P, a in sorted_zipped_list:
        l = w.length()
        aux = l - 2 * P.degree()
        row = [w.to_permutation(), w.reduced_word(), P, a,  aux,l]
        out.append(row)
        if aux==a: 
            D.append(w) 
            row = [f"\033[1;31m{cell}\033[0m" for cell in row]  # Red color
        table_data.append(row)

    headers = ['w', 'rex' , 'P(e,w)', 'a(w)', 'l(w) - 2 deg(P)', 'l(w)']

    print(f'\n\n\ntable of a-function values for S_{n}')
    print(tabulate(table_data, headers=headers, tablefmt='grid'))

    return out,D

data, D = sorted_result = a_function_table(n)






#######################################################################################

## --- filtration of the hecke algebra defined by KL basis elements having a-function values >= given value

#######################################################################################


class hecke_subquotient_via_lusztig_afunc():
    @staticmethod #### just to make the class independent of global functions
    def a_function(w):
        w = w.to_permutation()
        rsk=RSK(w)
        Q=rsk[1].conjugate()
        return sum([len(z)*(len(z)-1)//2 for z in Q])

    def __init__(self, W ,a_value):
        Welements = list(sorted(W.list(),key=lambda x: (x.length(),x.reduced_word())))
        self.weyl_group = W
        self.a_value = a_value
        self.W_elements = [w for w in Welements if hecke_subquotient_via_lusztig_afunc.a_function(w)==a_value]
        self.indices = [Welements.index(w) for w in self.W_elements]

        self.hecke_alg = IwahoriHeckeAlgebra(W, -scalar^2,1)
        self.KL_basis = self.hecke_alg.C()
        self.KL_basis_list = list(self.KL_basis.basis())

        self.dim = len(self.indices)

    def __repr__(self):
        return f'subquotient of the Hecke algebra of type A{self.weyl_group.rank()} having lusztigs a function value {self.a_value})'


    def coset_of(self,t):##### assume t is in the hecke algebra
        t = self.KL_basis(t)
        coeffs = t.to_vector()
        result_vec =  [0]*len(self.KL_basis_list)
        for i in self.indices:
            result_vec[i] = coeffs[i]
        out = 0
        for i,coeff in enumerate(result_vec):
            out += coeff * self.KL_basis_list[i]

        return out



sgn = hecke_subquotient_via_lusztig_afunc(W,0)
V = hecke_subquotient_via_lusztig_afunc(W,3)


##### can at least compare dimension with the specht mods
SGA = SymmetricGroupAlgebra(QQ, n)
partition = None 
for p in data:
    if p[3] == 3:
        partition = RSK(p[0])[1].conjugate().shape()

        break

SM = SGA.specht_module(partition)
print(SM.dimension()^2 == V.dim)

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