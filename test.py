from random import randint
from sympy import ntt,intt
import itertools
from typing import List, Tuple

f = open("output.txt","w")

n_set = [1024]
q = 12289

#seq1=[311, 239, 630, 441, 635, 70, 77, 96, 388, 561, 403, 313, 341, 554, 622, 445, 255, 292, 184, 487, 149, 21, 538, 142, 640, 151, 179, 538, 337, 550, 257, 517]
#seq2=[551, 411, 596, 263, 541, 202, 221, 312, 290, 88, 592, 625, 144, 487, 364, 430, 350, 128, 372, 599, 465, 202, 48, 9, 344, 375, 431, 396, 240, 397, 6, 466]

for n in n_set:
    for j in range(100):
        seq1=[randint(0,q-1) for x in range(n)]
        seq2=[randint(0,q-1) for x in range(n)]
        

        transform = []
        transform1 = ntt(seq1, q)
        transform2 = ntt(seq2,q)

        #print ("NTT1 :", transform1)
        for i in range(n):
            transform.append((transform1[i]*transform2[i]))
        #print(transform)
        transform = intt(transform,q)
        
        #print(transform1)
        #print(transform2)
        #print ("NTT :", transform)
        seq1=str(seq1)[1:-1]
        seq2=str(seq2)[1:-1]
        transform=str(transform)[1:-1]
        seq1=seq1.replace(',',"")
        seq2=seq2.replace(',',"")
        transform=transform.replace(',',"")
        f.write(str(seq1))
        f.write("\n")
        f.write(str(seq2))
        f.write("\n")
        f.write(str(transform))
        f.write("\n")
f.close()
