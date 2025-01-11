print("using RSK correspondence to calculate Lusztig's a function")

from collections import Counter

from tabulate import tabulate

R.<q> = LaurentPolynomialRing(QQ)


def sn_elements_one_line(n):
    return Permutations(n).list()

def n_choose2(n):
    return n*(n-1)//2

def a_function(w:list):
    rsk=RSK(w)
    Q=rsk[1].conjugate()
    return sum([n_choose2(len(z))for z in Q])


def get_multiplicities(lst):
    return Counter(lst)

def a_function_all(n):
        
    W = WeylGroup(f"A{n-1}")
    # [s1,s2,s3] = W.simple_reflections()
    KL = KazhdanLusztigPolynomial(W,q)
    elements= W.list()
    KLpolys = [KL.P(elements[0],s) for s in elements]
    sn_elements = [w.to_permutation() for w in elements]


    out= [a_function(w) for w in sn_elements]
    lengths = [Permutation(list(w)).number_of_inversions() for w in sn_elements]
    count = get_multiplicities(out)
    print('values and their multiplicities:')
    print(count)
    print('total number of values:', len(count))
    print(f'number of partitions of {n}:',number_of_partitions(n) )
    zipped_list = list(zip(sn_elements, KLpolys, out, lengths))
    sorted_zipped_list = sorted(zipped_list, key=lambda x: (x[-1], x[-2]))
    
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
    print(tabulate(table_data, headers=headers, tablefmt='grid'))

    return sorted_zipped_list,D




# Example usage
dara, D = sorted_result = a_function_all(3)